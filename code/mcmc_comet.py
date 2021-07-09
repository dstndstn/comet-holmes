# Copyright 2011 Dustin Lang (Princeton) and David W. Hogg (NYU).
# All rights reserved.

# BUGS:
# -----
# - no known bugs (well, search the code!)

import matplotlib
if __name__ == '__main__':
	matplotlib.use('Agg')
import multiprocessing
from glob import glob
import os
import datetime
import time
import sys
from math import pi
import numpy as np
from numpy.linalg import lstsq
import pylab as plt
if not hasattr(plt, 'tick_params'):
	print 'tick_params() is not available in this matplotlib version'
	def tick_params(**kwargs):
		pass
	plt.tick_params = tick_params
from matplotlib.patches import Ellipse, Polygon, Circle
from scipy import interpolate
from scipy.special import gammaln
#import markovpy
import emcee
from astrometry.util.file import *
from astrometry.util.multiproc import *
from astrometry.util import jpl
from astrometry.util import EXIF
from astrometry.util.ngc2000 import get_ngc
from astrometry.util.util import Tan
from astrometry.util.starutil_numpy import jdtomjd, radectoxyz, deg2distsq, datetomjd, mjdtodate, arcsec2distsq, xyztoradec, arcsec_between, ecliptic_basis
import astrometry.util.celestial_mechanics as cm

# in [AU^3 / yr^2]
GM_sun = cm.GM_sun

days_per_yr = 365.242199

# prior parameters
PRIOR_RADIUS = 1. # AU
PRIOR_ALPHA = 1.
PRIOR_BETA = 3.

def lnbeta(x, alpha, beta):
	return np.exp( gammaln(alpha + beta) - gammaln(alpha) - gammaln(beta)
		       + np.log(x) * (alpha-1.) + np.log(1.-x) * (beta-1.) )

def cosdeg(x):
	return np.cos(np.deg2rad(x))

# Returns mjd, E = (a,e,i,Omega,pomega,M)  in AU, radians
def EMB_ephem():
	jd,E = jpl.parse_orbital_elements(
		'''2454101.500000000 = A.D. 2007-Jan-01 00:00:00.0000 (CT)
 EC= 1.670361927937051E-02 QR= 9.832911829245575E-01 IN= 9.028642170169823E-04
 OM= 1.762399911457168E+02 W = 2.867172565215373E+02 Tp=  2454104.323526433203
 N = 9.856169820212966E-01 MA= 3.572170843984495E+02 TA= 3.571221738148128E+02
 A = 9.999947139070439E-01 AD= 1.016698244889530E+00 PR= 3.652534468934519E+02''',
									  needSystemGM=False)
	return jdtomjd(jd[0]),E[0]

# times: [days]
def EMB_xyz_at_times(times):
	t0,E = EMB_ephem()
	(a,e,I,Omega,pomega,M0,nil) = E
	# [radians / yr]
	dMdt = np.sqrt(GM_sun / a**3)
	# print 'Earth dMdt =', dMdt
	# ~= 2*pi
	obsxyz = []
	for t in times:
		M = M0 + dMdt * (t - t0) / days_per_yr
		(x,v) = cm.phase_space_coordinates_from_orbital_elements(
			a,e,I,Omega,pomega,M,GM_sun)
		obsxyz.append(x)
	obsxyz = np.array(obsxyz)
	return obsxyz

def inbox(x, y, xlo, xhi, ylo, yhi):
	return np.logical_and(np.logical_and(x >= xlo, x <= xhi),
						  np.logical_and(y >= ylo, y <= yhi))


class CometMCMCxv:
	#			  AU,  AU/year
	paramnames = ['x', 'v']
	def __init__(self):
		self.x = np.zeros(3)
		self.v = np.zeros(3)
		# epoch of the orbital elements (ie, M)
		# copy-n-pasted from JPL below
		self.epoch = jdtomjd(2454418.5)
		self.times = None
		self.earthxyz = None
		# spline subsampling
		self.Nspline = 200
		# self.pgood = 0.85
		# self.pexif = 0.75
		# self.centering = 2.5
		self.pgood = 0.8
		self.pexif = 0.8
		self.centering = 2.5
		self.set_times(datetomjd(datetime.datetime(2007, 7, 1)),
					   datetomjd(datetime.datetime(2008, 5, 1)),
					   0.05)

	def set_times(self, tmin, tmax, dtdays=None):
		self.tmax = tmax
		self.tmin = tmin
		if dtdays is None:
			dtdays = self.dtdays
			if dtdays is None:
				# default
				dtdays = 1.
		self.dtdays = dtdays
		self.times = np.arange(tmin, tmax, self.dtdays)
		self.set_ptime(np.ones_like(self.times))
		self.earthxyz = EMB_xyz_at_times(self.times)

	def set_ptime(self, pt):
		assert(pt.shape == self.times.shape)
		assert(np.sum(pt) > 0)
		self.ptime = pt.astype(float)
		self.ptime /= (self.ptime.sum() * self.dtdays)

	def set_ptime_from_samples(self, times, binwidth):
		bins = np.arange(self.tmin, self.tmax + binwidth, binwidth)
		H,xe = np.histogram(times, bins=bins)
		# np.digitize is 1-indexed
		inds = np.clip(np.digitize(self.times, bins) - 1, 0, len(bins)-1)
		# MAGIC +1: regularization
		self.set_ptime(H[inds].astype(float) + 1.)

	# [AU]
	def get_x(self):
		return self.x.copy()

	# [AU/yr]
	def get_v(self):
		return self.v.copy()

	# [AU]
	def set_x(self, x):
		self.x[:] = x[:]

	# [AU/yr]
	def set_v(self, v):
		self.v[:] = v[:]

	def set_epoch(self, mjd):
		self.epoch = mjd

	def set_params(self, p):
		self.x[:] = p[:3]
		self.v[:] = p[3:6]
		self.pgood = p[6]
		self.pexif = p[7]
		self.centering = np.exp(p[8])

	def get_params(self):
		return np.hstack((self.get_x(), self.get_v(),
				  self.pgood, self.pexif, np.log(self.centering)))

	def get_orbital_elements(self):
		return cm.orbital_elements_from_phase_space_coordinates(
			self.get_x(), self.get_v(), GM_sun)

	def set_params_from_jpl_string(self, s):
		# [AU], [AU/day], [jd]
		x,v,jd = jpl.parse_phase_space(s)
		x = x[0]
		# [AU/day] -> [AU/yr]
		v = v[0] * days_per_yr
		jd = jd[0]
		self.set_x(x)
		self.set_v(v)
		self.set_epoch(jdtomjd(jd))

	# year ~2000 visit
	def set_last_pass_params(self):
		s = '''2455638.500000000 = A.D. 2011-Mar-18 00:00:00.0000 (CT)
		-5.048281449175005E+00  4.365257402143395E-02 -9.429924038754340E-01
		9.663748139985810E-04 -5.542825718109530E-03 -1.422501670129140E-03
	    2.966181942514266E-02  5.135784829124947E+00 -7.358334729180353E-04
		'''
		self.set_params_from_jpl_string(s)
		
	def set_true_params(self):
		# from test-holmes-1.txt
		s = '''2454418.500000000 = A.D. 2007-Nov-14 00:00:00.0000 (CT)
		1.248613009072901E+00  2.025389080777020E+00  8.242272661173784E-01
		-7.973747689948190E-03  9.391385050634651E-03  1.214741189900433E-03
		1.454305558379576E-02  2.518052017168866E+00  3.997656275386785E-03
		'''
		self.set_params_from_jpl_string(s)

	def lnposterior(self, data):
		return self.lnprior() + self.lnlikelihood(data)

	# the prior is of the form
	# p(x,v,pgood,pexif) = p(x) p(v|x) p(pgood) p(pexif)
	# p(x) is a gaussian
	# p(v|x) is a beta distribution in v^2 between 0 and the unbinding v^2
	# p(pgood) = p(pexif) = uniform(0,1)
	def lnprior(self):
		lnp = 0.
		if self.tmax <= self.tmin:
			print 'lnprior(): tmax <= tmin -- punishing prior'
			lnp += -1000.
		if self.pgood <= 0 or self.pgood >= 1:
			print 'lnprior(): pgood=', self.pgood, '-- punishing prior'
			lnp += -1000.
		if self.pexif <= 0 or self.pexif >= 1:
			print 'lnprior(): pexif=', self.pexif, '-- punishing prior'
			lnp += -1000.
		if self.centering < 1.:
			print 'lnprior(): centering=', self.centering, '-- punishing prior'
			lnp += -1000.
		if self.get_energy() > 0:
			print 'lnprior(): unbound orbit -- punishing'
			lnp += -1000.
			return lnp # this must be here or else death and destruction
		x = self.get_x()
		lnp += -0.5 * np.dot(x,x) / PRIOR_RADIUS**2
		v2max = -2.0 * cm.potential_energy_from_position(x, GM_sun)
		v = self.get_v()
		lnp += lnbeta(np.dot(v,v) / v2max, PRIOR_ALPHA, PRIOR_BETA)
		# Jacobian factor converts [d^3v] to [dv^2]
		lnp += -1. * np.log(2. * pi * cm.norm(v))
		return lnp

	def get_energy(self):
		return cm.energy_from_phase_space_coordinates(
			self.get_x(), self.get_v(), GM_sun)

	def get_elements_at_time(self, t):
		C = list(self.get_orbital_elements()) + [GM_sun]
		# [cycles / yr]
		dMdt = np.sqrt(GM_sun / C[0]**3)
		#print 'comet dMdt:', dMdt
		# dMdt = 0.913 * 6.88 (yrs) = 2*pi (radians)
		C[5] += dMdt * ((t - self.epoch) / days_per_yr)
		return C

	# ecliptic coords, or equatorial?
	def heliocentric_xyz_at_times(self, times=None,
								  ecliptic=True):
		C = list(self.get_orbital_elements()) + [GM_sun]
		M0 = C[5]
		# [cycles / yr]
		dMdt = np.sqrt(GM_sun / C[0]**3)
		if times is None:
			times = self.times
		N = len(times)
		XYZ = np.empty((N,3))
		sun = np.zeros(3)
		for i,t in enumerate(times):
			C[5] = M0 + dMdt * ((t - self.epoch) / days_per_yr)
			if ecliptic:
				XYZ[i,:] = cm.orbital_elements_to_xyz(C, sun,
													  light_travel=False,
													  normalize=False)
			else:
				(x,dx) = cm.orbital_elements_to_ss_xyz(C, sun,
													   light_travel=False)
				XYZ[i,:] = dx
		return XYZ
		

	# Note, this returns "delta-xyz" from EMB, in Equatorial coords.
	# times are mjds
	def xyz_at_times(self, times=None, light_travel=True):
		C = list(self.get_orbital_elements()) + [GM_sun]
		M0 = C[5]
		# [cycles / yr]
		dMdt = np.sqrt(GM_sun / C[0]**3)
		if times is None:
			times = self.times
			earthxyz = self.earthxyz
			Nsp = self.Nspline
		else:
			Nsp = 1
			earthxyz = EMB_xyz_at_times(times)
		while len(times)/Nsp < 6 and Nsp > 1:
			Nsp = max(1, int(Nsp / 2))

		SS = slice(0, len(times), Nsp)
		N = len(times[SS])

		XYZ = np.empty((N,3))
		for i,(ex,t) in enumerate(zip(earthxyz[SS], times[SS])):
			C[5] = M0 + dMdt * ((t - self.epoch) / days_per_yr)
			XYZ[i,:] = cm.orbital_elements_to_xyz(C, ex, light_travel)
		
		if Nsp > 1:
			(S,u) = interpolate.splprep([XYZ[:,0], XYZ[:,1], XYZ[:,2]], u=times[SS], s=0)
			X,Y,Z = interpolate.splev(times, S)
			XYZ = np.vstack((X,Y,Z)).T
			
		# just in case any of our shih is hucked up.
		XYZ /= np.sqrt(np.sum(XYZ**2, axis=1))[:,np.newaxis]
		return XYZ

	def radec_at_times(self, times=None, light_travel=True):
		xyz = self.xyz_at_times(times, light_travel)
		return xyztoradec(xyz)

	def find_points_in_wcs(self, wcs, xyz,
						   xlo=None, xhi=None, ylo=None, yhi=None):
		'''
		Returns the set of indices in "xyz" that are inside the
		given WCS, and optionally inside the given ROI.

		xyz has shape (N,3)
		'''
		if xlo is None:
			xlo = 0.5
		if xhi is None:
			xhi = wcs.imagew + 0.5
		if ylo is None:
			ylo = 0.5
		if yhi is None:
			yhi = wcs.imageh + 0.5

		# Pre-filter based on WCS center and radius
		xc,yc = (xlo+xhi)/2., (ylo+yhi)/2.
		xyzc = np.array(wcs.pixelxy2xyz(xc, yc))
		r2 = arcsec2distsq(wcs.pixel_scale() *
						   np.hypot(xhi-xlo, yhi-ylo)/2.)
		R2 = ((xyz - xyzc[np.newaxis,:])**2).sum(axis=1)
		# make the pre-filter conservative by fudge factor
		fudge = 1.05
		I = np.flatnonzero(R2 < (r2 * fudge))
		if len(I) == 0:
			return np.array([])

		# Push the points that passed the pre-filter through the WCS
		x,y = wcs.xyz2pixelxy(xyz[I,:])
		x,y = np.array(x), np.array(y)
		J = inbox(x, y, xlo, xhi, ylo, yhi)
		return I[J]

	# God's pwn function
	# MAGIC: 0.5 days = 12 hrs = span of time zones
	def find_times_slice_within_12_hours_of_mjd(self, time):
		return (self.times > (time - 0.5)) * (self.times < (time + 0.5))

	# NOTE: this likelihood assumes that self.times are evenly spaced.
	def lnlikelihood(self, data, getfgbgs=False):
		try:
			cxyz = self.xyz_at_times()
		except:
			return -100. * len(data)
		pgood = np.clip(self.pgood, 0.001, 0.999)
		pexif = np.clip(self.pexif, 0.001, 0.999)
		# MAGIC 1/(4*pi..yadda) = the sky in arcsec^2
		bg = (1. - pgood) / (4. * pi * 206265.**2)
		if self.centering > 1:
			sfc = np.sqrt(float(self.centering))
			f1 = 0.5 - 0.5 / sfc
			f2 = 0.5 + 0.5 / sfc
		else:
			f1 = 0.
			f2 = 1.
		lnp = 0.
		wcs = Tan()
		if getfgbgs:
			fgbgs = []
		for i, (wcsvals, date) in enumerate(data):

			# Compute time prior p(time)
			# MAGIC: we set bad dates to zero
			if date > 1:
				H = self.find_times_slice_within_12_hours_of_mjd(date)
				if H.sum() > 0:
					ptime = (1. - pexif) * self.ptime
					# BUG / MAGIC: the following lines rely on dtdays being in days!
					ptime[H] += pexif
					ptime /= (np.sum(ptime) * self.dtdays)
			else:
				ptime = self.ptime

			# Count time steps that lie inside this image's WCS (with ln(centering))
			wcs.set(*wcsvals)
			xlo,xhi = f1 * wcs.imagew + 0.5, f2 * wcs.imagew + 0.5
			ylo,yhi = f1 * wcs.imageh + 0.5, f2 * wcs.imageh + 0.5
			J = self.find_points_in_wcs(wcs, cxyz, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi)

			# Could:
			# look for places where J goes from 0 to 1 (ingress) and 1 to 0 (egress)
			# find crossing point on image, convert to time, and scale this timestep
			# (binary search for Xgress time)

			# --not quite -- J is a list of indices -- so have to extend out from non-contiguous endpoints.

			# Integrate p(time)
			if len(J) == 0:
				totalp = 0.
			else:
				totalp = ptime[J].sum()

				(N,three) = cxyz.shape
				# indices of left edge of breaks (ingresses & egresses)
				breaks = []
				if J[0] > 0:
					breaks.append(J[0]-1)
				if J[-1] < (N-1):
					breaks.append(J[-1])
				jumps = np.flatnonzero(np.diff(J) - 1)
				for j in jumps:
					breaks.append(J[j])
					breaks.append(J[j+1]-1)

				#print 'breaks:', breaks
				for jlo in breaks:
					xyzlo = cxyz[jlo,:]
					xyzhi = cxyz[jlo+1,:]
					dxyz = xyzhi - xyzlo
					x,y = wcs.xyz2pixelxy(xyzlo)
					inlo = inbox(x, y, xlo, xhi, ylo, yhi)
					x,y = wcs.xyz2pixelxy(xyzhi)
					inhi = inbox(x, y, xlo, xhi, ylo, yhi)
					assert(inlo or inhi)
					assert(not(inlo and inhi))
					assert(inlo != inhi)
					# binary search
					flo, fhi = 0., 1.
					for step in range(10):
						fmid = (flo + fhi)/2.
						#print 'flo fmid fhi', flo, fmid, fhi
						if flo >= fhi:
							break
						xyz = xyzlo + fmid * dxyz
						x,y = wcs.xyz2pixelxy(xyz)
						inmid = inbox(x, y, xlo, xhi, ylo, yhi)
						if inmid == inlo:
							# move up lo
							flo = fmid
						else:
							# move down hi
							fhi = fmid
					#print 'fmid', fmid

					## ISSUE: we want to integrate ptime over time.
					# Our samples are *midpoints* in ptime histogram bins.
					# If a sample clips the edge of an image, then it
					# gets overcounted (half its ptime bin is outside the
					# image).

					# Since samples are midpoints of ptime bins, the bin we
					# integrate from changes at fmid=0.5.
					
					if inlo:
						# This is an egress: jlo is inside
						# (and has already been counted).
						if fmid < 0.5:
							# part of ptime bin "jlo" was counted but
							# it's actually outside the image (fmid to 0.5)
							totalp -= ptime[jlo] * (0.5 - fmid)
						else:
							# part of the next ptime bin should have been
							# counted.
							totalp += ptime[jlo+1] * (fmid - 0.5)
					else:
						# This is an ingress: jlo is outside and jhi=jlo+1 is inside
						# (and has already been counted)
						if fmid < 0.5:
							# part of the jlo bin should be counted (from fmid to 0.5).
							totalp += ptime[jlo] * (0.5 - fmid)
						else:
							# part of jhi bin was overcounted (from 0.5 to fmid)
							totalp -= ptime[jlo+1] * (fmid - 0.5)

			fg = totalp * self.dtdays * pgood / (wcs.pixel_scale()**2 * (yhi-ylo)*(xhi-xlo))
			# print 'lnlikelihood(): sum(J), 'in image, totalp=', totalp, 'fg', fg, 'bg', bg
			if getfgbgs:
				fgbgs.append((fg,bg))
			lnp += np.log(fg + bg)
		if getfgbgs:
			return lnp,fgbgs
		return lnp


# a function-like object that encapsulates the data we need to do
# a lnposterior call and can be passed to multiprocessing.
class CometLnPosterior(object):
	def __init__(self, C, data):
		self.C = C
		self.data = data
	def __call__(self, params):
		#print self.C
		self.C.set_params(params)
		return self.C.lnposterior(self.data)


def plot_traj(spargs, pos, data2,
			  epoch, fn, fn2, color):
	print 'Plotting trajectory', fn
	plt.clf()
	plt.subplots_adjust(**spargs)
	plot_trajectory_2(pos, data2, epoch=epoch, color=color)
	# Check -- plot the JPL ephemeris in RA,Dec
	# this is from Earth (399) not EMB
	# Looks good!
	#(ra,dec,jd) = jpl.parse_radec(open('holmes-ephem.txt').read())
	#plt.plot(ra, dec, 'b-', alpha=0.5)
	#plt.axis([70, 35, 30, 60])
	plt.savefig(fn)
	if fn2:
		plt.savefig(fn2)
	print 'Done plotting trajectory', fn


# For multiprocessing -- called by each worker in the Pool when it
# starts.
def global_init():
	import matplotlib
	matplotlib.rc('font', family='computer modern roman')
	matplotlib.rc('font', size=14)
	matplotlib.rc('text', usetex=True)
	matplotlib.rc('text.latex', preamble=r'\usepackage{bm}\usepackage{amsmath}')
	#dpi = 300
	dpi = 100
	# figure(dpi=) does nothing.  Use savefig(dpi=) instead.
	matplotlib.rc('savefig', dpi=dpi)
	plt.figure(figsize=(6,6)) #, dpi=dpi)
	


# 2012-02-02:
#    mkdir 2012-02-02
#    mv mcmc-default-nw064-bw008-emp-*.pickle 2012-02-02/
#    mv traj-default-nw064-bw008-emp-*.png 2012-02-02/
#    python -u mcmc_comet.py --maxiter=5001 --emp-r=1 > log 2>&1 &
# I stopped this at step ~2874 in order to simplify the initialization steps
# (to make the paper cleaner); moved results to:
# 2012-02-03 .  That run was started around rev 20230 (bit earlier, but no
# substantive changes mcmc_comet.py anyway);
#    mkdir 2012-02-03
#    mv mcmc*.pickle 2012-02-03/
#    mv traj-default-nw064-bw008-*.png 2012-02-03/
#    mv log 2012-02-03/
#
# At rev 20231:
#    nice -n 20 python -u mcmc_comet.py --maxiter=5001 --emp-r=1 >> log 2>&1 &
#
# python -u mcmc_comet.py --maxiter=5001 --emp-r=1 >> log 2>&1 &
# 
#

#
# [older:]
# paper:
#
#   # python mcmc_comet.py --maxiter=2001 --emp-r=1
# (no, actually just:)
#   python mcmc_comet.py --maxiter=2001
#
#   python mcmc_comet.py --maxiter=2001 --exif-plots
#
# and for checking....
#   python mcmc_comet.py --maxiter=2001 --prob-plots
#
def main():
	import optparse
	parser = optparse.OptionParser()
	parser.add_option('--small-only', dest='smallonly', default=False, action='store_true', help='Use only (angularly) small images')
	parser.add_option('--hack-out-crap', '--hoc', dest='hoc', default=False, action='store_true', help='Remove images that might be screwing us')
	parser.add_option('--threads', dest='threads', default=16, type=int, help='Use this many concurrent processors')
	parser.add_option('--walkers', dest='walkers', default=64, type=int, help='Use this many MCMC walkers')
	parser.add_option('--binwidth', dest='binwidth', default=8, type=int, help='Make cheater prior with this bin width (days)')
	parser.add_option('--maxiter', dest='maxiter', default=500, type=int, help='Do no more than this number of iterations')
	parser.add_option('--exif-plots', dest='exifplots', default=False, action='store_true', help='Make EXIF plots and exit.')
	parser.add_option('--emp-r', dest='emprad', type=float, default=None, help='Radius for empirical initialization')
	parser.add_option('--prob-plots', dest='probplots', default=False, action='store_true', help='Make probability-surface plots and exit.')
	opt,args = parser.parse_args()
	np.random.seed(42)

	print 'main(): Munging data...'
	print '  reading WCS...'
	wcsfns = glob('2010-04-16-holmes-dedup/holmes-*.wcs')
	wcsfns.sort()
	wcs = [Tan(fn,0) for fn in wcsfns]

	# ad-hockery!
	if opt.hoc:
		keep = []
		for i,w in enumerate(wcs):
			rmin,rmax,decmin,decmax = w.radec_bounds()
			#print 'RA', rmin, rmax, 'Dec', decmin, decmax
			if rmin > 40. and rmax < 44. and decmin > 30. and decmax < 35.:
				print 'main(): removing image', wcsfns[i].replace('.wcs','.jpg')
			else:
				keep.append(i)
		wcsfns = [wcsfns[i] for i in keep]
		wcs = [wcs[i] for i in keep]

	# if asked, keep only (angularly) smallest half
	if opt.smallonly:
		scales = np.array([w.pixel_scale()*np.sqrt(w.imagew*w.imageh) for w in wcs])
		I = np.argsort(scales)[:len(wcs)/2]
		wcs = [wcs[i] for i in I]

	# parse dates from EXIF data in jpegs
	dates,exifobjs = parse_exif([fn.replace('.wcs','.jpg')
								 for fn in wcsfns])

	gooddates = [d for d in dates if d is not None]
	print 'main(): Got', len(gooddates), 'dates from EXIF'
	goodmjds = [datetomjd(d) for d in gooddates]

	data = []
	for w,d in zip(wcs, dates):
		mjd = 0
		if d is not None:
			mjd = datetomjd(d)
		data.append(((w.crval[0], w.crval[1], w.crpix[0], w.crpix[1],
					  w.cd[0], w.cd[1], w.cd[2], w.cd[3],
					  w.imagew, w.imageh),
					 mjd))
	print 'main(): burp.'
	print 'main(): Number of images:', len(data)

	# make comet object
	C = CometMCMCxv()
	C.set_ptime_from_samples(goodmjds, opt.binwidth)
	# we duplicate this comet object a lot of times by using
	# multiprocessing -- what bits of data do we actually need
	# to share? -- (These are all constant)
	#  .epoch, .times, .earthxyz, .Nspline
	#  .tmax, .tmin, .dtdays, .ptime
	# Could use multiprocessing.Value/Array for these.
	lnprobfunc = CometLnPosterior(C, data)

	if opt.exifplots:
		from exifplots import exifplots
		C.set_true_params()
		exifplots(C, data, exifobjs)
		sys.exit(0)

	epoch = None
	if opt.emprad:
		epoch,cx0,cv0 = empirical_init(opt.emprad, goodmjds, dates, data, wcs)
		C.set_epoch(epoch)
		C.set_x(cx0)
		C.set_v(cv0)

		print 'Epoch', epoch
		
	else:
		C.set_true_params()

	# MAGIC step sizes
	#    xxx, vvv, pgood, pexif, eta
	#stepsizes = np.array([1e-5]*3 + [5e-4]*3 + [0.001, 0.003, 0.001])
	stepsizes = np.array([1e-5]*3 + [1e-3]*3 + [0.01, 0.01, 0.01])

	mp = multiproc(opt.threads, global_init, ())
	if opt.threads == 1:
		global_init()

	if opt.probplots:
		p0 = C.get_params().copy()
		for i,(step,ll) in enumerate(zip(stepsizes, ['x0 [AU]','x1 [AU]','x2 [AU]',
													 'v0 [AU/yr]','v1 [AU/yr]','v2 [AU/yr]',
													 'pgood','pexif','centering'])):
			print 'Making ln-posterior plot for parameter', i
			pvals = np.linspace(p0[i] - 3.*step, p0[i] + 3.*step, 31)
			p = p0.copy()
			allparams = []
			for pval in pvals:
				p[i] = pval
				allparams.append(p.copy())
			lnps = mp.map(lnprobfunc, allparams)

			plt.clf()
			plt.plot(pvals, lnps - lnps[len(lnps)/2], 'r-')
			plt.xlabel(ll)
			plt.ylabel('delta ln(posterior)')
			plt.axvline(p0[i], color='k', alpha=0.5)
			plt.axvline(p0[i]-step, color='k', alpha=0.5)
			plt.axvline(p0[i]+step, color='k', alpha=0.5)
			plt.xlim(pvals[0], pvals[-1])
			plt.savefig('param%i.png' % i)

		sys.exit(0)

	# What are we going to name our output files?
	suffix = ''
	if opt.smallonly: suffix += '-so'
	if opt.hoc: suffix += '-hoc'
	if suffix == '': suffix += '-default'
	suffix += '-nw%03d' % opt.walkers
	suffix += '-bw%03d' % opt.binwidth
	if opt.emprad == 1.:
		suffix += '-emp'
	elif opt.emprad is not None:
		suffix += '-emp%.2f' % opt.emprad


	# Initialize MCMC sampler...
	nwalkers = opt.walkers
	ndim	 = len(C.get_params())
	# randomize walker positions around the initialization point.
	initpos = []
	for i in range(nwalkers):
		p = C.get_params()
		p += np.random.normal(size=p.shape) * stepsizes
		initpos.append(p)

	sampler = markovpy.EnsembleSampler(nwalkers, ndim, lnprobfunc,
									   pool=mp.pool)

	# Get true parameters to set plot limits...
	Ctrue = CometMCMCxv()
	Ctrue.set_true_params()
	ptrue = Ctrue.get_params()
	del Ctrue
	postrue = ptrue[:6]
	scale = np.array([1,1,1,365.25,365.25,365.25])
	lims = dict([(i,[postrue[i]*scale[i] - 3e-1,
					 postrue[i]*scale[i] + 3e-1])
				  for i in range(len(scale))])
	lims2 = dict([(0,[0., 1.]), (1, [0., 1.]), (2, [0., 2.])])

	labels2 = ['pgood', 'pexif', 'ln(centering)']
	#spargs = dict(bottom=0.15, left=0.15, right=0.95, top=0.95)
	spargs = dict(bottom=0.1, left=0.1, right=0.98, top=0.98)
	
	# Initialize MCMC sampler state:
	lnprob,state = None,None
	pos = np.array(initpos)

	for k in range(opt.maxiter):
		if (k%10 == 0):
			for color in [True,False]:
				if color:
					bwstr = ''
				else:
					bwstr = '-bw'
				fn = 'traj%s-%04d%s.png' % (suffix, k, bwstr)
				if os.path.exists(fn):
					print 'Plot file already exists:', fn
				else:
					fn2 = None
					if (k%100 == 0):
						fn2 = 'traj%s-%04d%s.pdf' % (suffix, k, bwstr)
					mp.apply(plot_traj, (spargs, pos, data,
										 epoch, fn, fn2, color))

		if False:
			fn = 'manyd%s-%04d.png' % (suffix, k)
			if (k%10 == 0) and (not os.path.exists(fn)):
				print 'Plotting distribution', fn
				x = mp.apply(applyplot,
							 (spargs, 6, pos[:,:6]*scale, ptrue[:6]*scale, lims, None, fn))

			fn = 'manyd2%s-%04d.png' % (suffix, k)
			if (k%10 == 0) and (not os.path.exists(fn)):
				print 'Plotting distribution', fn
				x = mp.apply(applyplot,
							 (spargs, 3, pos[:,6:9], ptrue[6:9], lims2, labels2, fn))

		fn0 = 'mcmc%s-%03i.pickle' % (suffix, k)
		fn = 'mcmc%s-%04i.pickle' % (suffix, k)
		if os.path.exists(fn):
			print 'Reading pickle', fn
			(nil,pos,lnprob,state) = unpickle_from_file(fn)
		elif os.path.exists(fn0):
			print 'Reading pickle', fn0
			(nil,pos,lnprob,state) = unpickle_from_file(fn0)
		else:
			print 'Running MCMC', fn
			t0 = datetime.datetime.now()
			niter = 1
			pos,lnprob,state = sampler.run_mcmc(pos, state, niter, lnprobinit=lnprob)
			print 'main(): current acceptance fraction:', 100. * sampler.acceptance_fraction
			t1 = datetime.datetime.now()
			dt = t1-t0
			dt = (dt.microseconds + (dt.seconds + dt.days * 24 * 3600.) * 1e6) / 1e6
			print 'MCMC took', dt, 'sec'
			pickle_to_file((k,pos,lnprob,state), fn)

	mp.waitforall()


def plot_trajectory_2(pos, data, npts=16,
					  S=[], epoch=None, color=True):
	I = plt.imread('holmes-footprint-1.png')
	plotwcs = Tan('holmes-footprints-1.wcs',0)
	H,W,four = I.shape

	# dim it a bit?
	#I[:,:,:3] *= 0.95
	I[:,:,:3] *= 0.8

	# Flipped vertically...
	plt.imshow(I, interpolation='nearest',
			   # FITS WCS coords
			   extent=(1,W,H,1))
	ax = plt.axis()
	
	# Plot sample trajectories
	c = CometMCMCxv()
	if epoch is not None:
		c.set_epoch(epoch)
	c.dtdays = 1
	c.Nspline = 10

	# linecolor
	if color:
		#lc = (0,1,0)
		lc = (0.1,1.0,0.1)
		# lc = (0,0.9,0)
		#lc = (0,0.3,0)
	else:
		lc = (1,1,1)
	
	p1 = None
	samplealpha = 0.5
	for posi in pos[:npts]:
		c.set_params(posi)
		(ras,decs) = c.radec_at_times()
		x,y = plotwcs.radec2pixelxy(ras, decs)
		p1 = plt.plot(x, y, '-', color=lc, alpha=samplealpha, linewidth=0.5)
		# dots on endpoints
		plt.plot([x[0], x[-1]], [y[0], y[-1]], '.', color=lc, alpha=samplealpha)

		# Fake plot with stronger alpha for the legend.
		p1 = plt.plot([ax[1]+1, ax[1]+2], [0,0], '-', color=lc, alpha=1., linewidth=1)


	# Including epoch!
	c.set_true_params()
	(ras,decs) = c.radec_at_times()
	x,y = plotwcs.radec2pixelxy(ras, decs)
	p2 = plt.plot(x, y, ':', color=lc, alpha=1, linewidth=2.0)

	# Label some dates
	for day,ha,va in [ #((2007,  8, 1), 'left',  'bottom'),
		((2007,  8, 5), 'left',  'bottom'),
		((2007, 10, 1), 'right', 'top'),
		((2007, 12, 1), 'left',  'bottom'),
		((2008,  2, 1), 'left',  'top'),
		((2008,  3, 25), 'left',  'top') ]:
		#((2008,  4, 1), 'left',  'top') ]:
					   
		D = datetime.datetime(*day)
		d = datetomjd(D)
		ra,dec = c.radec_at_times(np.array([d]))
		x,y = plotwcs.radec2pixelxy(ra, dec)
		plt.plot(x, y, 'o', mec=lc, mfc=lc)#'none')
		# offset
		dr = 0.5 * (1. if ha == 'right' else -1.)
		dd = 0.5 * (1. if va == 'bottom' else -1.)
		x,y = plotwcs.radec2pixelxy(ra+dr, dec+dd)
		plt.text(x, y, str(D.date()), color=lc, size=10,
				 ha=ha, va=va, 
				 bbox=dict(facecolor='k', ec=lc, alpha=0.25))

	if p1 is not None:
		L = plt.legend((p2,p1), ('JPL', 'samples'), prop=dict(size=12))
		#if not color:
		if True:
			# Make the legend have a dark background
			F = L.get_frame()
			F.set_facecolor('0.2')

		# Set legend text = lc
		for t in L.get_texts():
			t.set_color(lc)

	plt.xlabel('RA (deg)')
	plt.ylabel('Dec (deg)')
	#plt.axis([70,35,31,59])

	# Mark the California nebula = NGC1499
	ngc = get_ngc(1499)
	x,y = plotwcs.radec2pixelxy(ngc.ra, ngc.dec)
	pixscale = plotwcs.pixel_scale()
	R = ngc.size * 60. / pixscale / 2.
	#plt.plot(x, y, 'o', mec=lc, mfc=lc)#'none')
	plt.gca().add_artist(Circle(xy=(x,y), radius=R, ec=lc, fc='none', alpha=1, lw=1))
	plt.text(x + R + 10, y, 'NGC 1499', color=lc, va='center', ha='left', fontsize=10)

	# DIY grid
	ras = np.linspace(30., 75., 500)
	decs = np.linspace(25., 65., 500)
	gc = (0,0.4,0)
	tv,tt = [],[]
	for ra in np.arange(30, 75, 5):
		R = ra + np.zeros_like(decs)
		x,y = plotwcs.radec2pixelxy(R, decs)
		plt.plot(x, y, color=gc, alpha=0.3)
		# Decs are in increasing order;
		# Bottom of the plot window is WCS y=H
		I = np.flatnonzero(y < H)[0]
		tt.append(ra)
		tv.append(x[I])
	plt.xticks(tv,tt)
	tv,tt = [],[]
	for dec in np.arange(30, 60, 5):
		D = dec + np.zeros_like(ras)
		x,y = plotwcs.radec2pixelxy(ras, D)
		plt.plot(x, y, color=gc, alpha=0.3)
		# RAs go right-to-left in the plot
		I = np.flatnonzero(x < 0)[0]
		tt.append(dec)
		tv.append(y[I])
	plt.yticks(tv,tt)

	# MAGIC 0.5s etc
	x0 = 20
	y0 = 990
	dy = -30
	tc = lc
	nhyp = 10
	hstyle = dict(color=tc, alpha=samplealpha, fontsize=10)
	I = np.arange(nhyp)

	pgood = pos[I,6]
	pexif = pos[I,7]
	eta = 1./np.exp(pos[I,8])

	for i,(sym,vals) in enumerate([(r'p_{\mathrm{good}}', pgood),
								  (r'p_{\mathrm{EXIF}}',  pexif),
								  (r'\eta',               eta),]):
		plt.text(x0, y0 + i*dy,
				 r'($\boldsymbol{%s}$ {\small %s})' % (sym, ' '.join(['%.2f' % v for v in vals])),
				 **hstyle)
	plt.axis(ax)

def applyplot(spargs, nd, pos, true, lims, labels, fn):
	plt.clf()
	plt.subplots_adjust(**spargs)
	manyd_plot(nd, pos, true, lims, labels=labels)
	plt.savefig(fn)
	print 'saved', fn

def manyd_plot(ndim,pos,ptrue,lims={},labels=None):
	if labels is None:
		labels = ['p%02d' % i for i in range(ndim)]
	for i in range(ndim):
		for j in range(ndim):
			plt.subplot(ndim, ndim, i*ndim + j + 1)
			plt.axvline(ptrue[j], color='g', alpha=0.5, lw=0.5)
			if i == j:
				plt.hist(pos[:,i], 20, range=lims.get(i, None))
			else:
				plt.axhline(ptrue[i], color='g', alpha=0.5, lw=0.5)
				plt.plot(pos[:,j], pos[:,i], 'k.', alpha=0.25)
			if j in lims:
				plt.xlim(*lims[j])
			else:
				lims[j] = plt.xlim()
			if i != j:
				if i in lims:
					plt.ylim(*lims[i])
				else:
					lims[i] = plt.ylim()
			plt.xlabel('')
			plt.ylabel('')
			plt.tick_params(labelbottom=False, labelleft=False,
							labelsize=8)
			if i == 0:
				plt.tick_params(labeltop=True)
			if i == ndim-1:
				plt.tick_params(labelbottom=True)
				plt.xlabel(labels[j])
			if j == 0:
				plt.tick_params(labelleft=True)
				plt.ylabel(labels[i])
			if j == ndim-1:
				plt.tick_params(labelright=True)
			if i == j:
				plt.tick_params(labelleft= False, labelright=False)
				plt.ylabel('')

def parse_exif(jpegfns):
	print '  reading EXIF...'
	dates = []
	exifobjs = []
	for imgfn in jpegfns:
		exif = EXIF.process_file(open(imgfn), details=False)
		exifobjs.append((imgfn, exif))
		imgdate = exif.get('EXIF DateTimeOriginal')
		if not imgdate:
			imgdate = exif.get('Image DateTime')
		if imgdate:
			dates.append(datetime.datetime.strptime(str(imgdate),
													'%Y:%m:%d %H:%M:%S'))
		else:
			dates.append(None)
	return dates,exifobjs

def empirical_init(rad, goodmjds, dates, data, WCS):
	print 'Initializing from the data with radius', rad
	# Initialize orbit using images close to the median date.
	medmjd = np.median(goodmjds)
	nearby = (np.abs(goodmjds - medmjd) < 7)
	print 'number of images within 7 days:', sum(nearby)
	# Move the epoch to the nearest integer MJD.
	epoch = medmjd
	epoch = float(round(epoch))
	print 'epoch', epoch

	# True radius
	# embxyz0 = EMB_xyz_at_times(np.array([epoch]))
	# embxyz0 = embxyz0[0,:]
	# print 'emb xyz0:', embxyz0
	# E = list(C.get_orbital_elements()) + [GM_sun]
	# (cxyz0,v0) = cm.phase_space_coordinates_from_orbital_elements(*E)
	# print 'comet xyz0:', cxyz0
	# print 'R:', np.sqrt(np.sum((embxyz0 - cxyz0)**2)), 'AU'

	Igooddates = np.array([i for i,d in enumerate(dates) if d is not None])
	Inear = Igooddates[nearby]
	nearwcs = [WCS[i] for i in Inear]
	nrd = np.array([wcs.radec_center() for wcs in nearwcs])
	nra  = nrd[:,0]
	ndec = nrd[:,1]
	# [deg]
	nrad = np.array([wcs.pixel_scale() * np.hypot(wcs.imagew, wcs.imageh)
					 for wcs in nearwcs]) / 3600.
	# cut on RA,Dec
	mra = np.median(nra)
	mdec = np.median(ndec)
	# MAGIC: search radius in deg
	Jnear = np.hypot(nra - mra, ndec - mdec) < 5.
	Inear = Inear[Jnear]
	nra = nra[Jnear]
	ndec = ndec[Jnear]
	nrad = nrad[Jnear]
	nt = np.array(goodmjds)[nearby][Jnear]

	# least-squares fit for dRA,dDec / dt
	N = len(nra)
	print 'Initializing with', N, 'nearby images'
	# weight ~ 1 / (image size)
	sqrtW = 1./nrad
	X = np.zeros((N,2))
	X[:,0] = 1. * sqrtW
	X[:,1] = (nt-epoch) * sqrtW
	y = np.zeros((N,2))
	y[:,0] = nra * sqrtW
	y[:,1] = ndec * sqrtW
	(b,resid,rank,eigs) = lstsq(X, y)
	print 'b', b
	ra0,dec0 = b[0,:]
	dradt,ddecdt = b[1,:]

	if False:
		plt.clf()
		plt.subplot(2,1,1)
		plt.errorbar(nt-epoch, nra, yerr=nrad, fmt=None)
		a = plt.axis()
		tt = np.array(a[:2])
		plt.plot(tt, ra0 + dradt * tt, 'b-')
		plt.xlim(a[0],a[1])
		plt.ylabel('RA (deg)')
		plt.subplot(2,1,2)
		plt.errorbar(nt-epoch, ndec, yerr=nrad, fmt=None)
		plt.plot(tt, dec0 + ddecdt * tt, 'b-')
		plt.xlim(a[0],a[1])
		plt.ylabel('Dec (deg)')
		plt.xlabel('dt (days)')
		plt.savefig('ravst.png')

	# Find velocity and compute orbital elements from x,v.
	dt = 1.
	ra1  = ra0  + dradt  * dt
	dec1 = dec0 + ddecdt * dt
	xyz0 = radectoxyz(ra0, dec0)[0]
	xyz1 = radectoxyz(ra1, dec1)[0]
	#print 'xyz:', xyz0, xyz1
	# xyz are in celestial coords; convert to solar-system coords
	(antieq, antisol, antipole) = ecliptic_basis(eclipticangle = -23.439281)
	xyz0 = xyz0[0] * antieq + xyz0[1] * antisol + xyz0[2] * antipole
	xyz1 = xyz1[0] * antieq + xyz1[1] * antisol + xyz1[2] * antipole
	# Find observer (EMB) positions at epoch
	embxyz = EMB_xyz_at_times(np.array([epoch, epoch + dt]))
	embx0 = embxyz[0,:]
	embx1 = embxyz[1,:]
	# Put comet at radius R from the EMB and find its x,v
	R = rad
	cx0 = embx0 + R * xyz0
	cx1 = embx1 + R * xyz1
	# AU/day
	cv0 = (cx1 - cx0) / dt
	# -> AU/yr
	cv0 *= days_per_yr
	return epoch, cx0, cv0

def eyeball_hyperparams():
	if False:
		for centering in np.arange(1.,4.01,0.25):
			C.set_ptime_from_samples(goodmjds, 8.)
			C.centering = centering
			print centering, C.lnposterior(data)
		sys.exit(0)

	if False:
		for binwidth in np.arange(5.,13.01):
			C.centering = 2.5
			C.set_ptime_from_samples(goodmjds, binwidth)
			print binwidth, C.lnposterior(data)
		sys.exit(0)

	if False:
		C.centering = 2.5
		C.set_ptime_from_samples(goodmjds, 8.)
		for pgood in np.arange(0.05,1.,0.1):
			C.pexif = 0.75
			C.pgood = pgood
			print C.pgood, C.pexif, C.lnposterior(data)
		for pexif in np.arange(0.05,1.,0.1):
			C.pgood = 0.85
			C.pexif = pexif
			print C.pgood, C.pexif, C.lnposterior(data)
		sys.exit(0)



if __name__ == '__main__':
	main()
