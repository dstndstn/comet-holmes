import matplotlib
matplotlib.use('Agg')

import cPickle as pickle
import numpy as np

from numpy import *
from numpy.linalg import *
from numpy.random import normal
from astrometry.util.mcmc import *
from astrometry.util.file import *
from astrometry.util.pyfits_utils import *
from astrometry.util.starutil_numpy import *
from astrometry.util import jpl
from matplotlib.patches import Circle,Ellipse
from celestial_mechanics import *
from comet_likelihood import *
import datetime
import sys
import numpy
import time
from pylab import *
import pylab

def savefig(x):
	print 'saving', x
	pylab.savefig(x)

def norm(x):
	return sqrt(dot(x, x))

def cosdeg(x):
	return cos(deg2rad(x))

# approx.
GM_sun = 2.95912e-04 # AU^3/d^2

def daterange(start, end, step):
	'''
	'''
	dt = datetime.timedelta(step)
	times = []
	t = start
	while t < end:
		times.append(t)
		t += dt
	return times

# Returns mjd,E
def EMB_ephem():
	jd,E = jpl.parse_orbital_elements('''2454101.500000000 = A.D. 2007-Jan-01 00:00:00.0000 (CT)
 EC= 1.670361927937051E-02 QR= 9.832911829245575E-01 IN= 9.028642170169823E-04
 OM= 1.762399911457168E+02 W = 2.867172565215373E+02 Tp=  2454104.323526433203
 N = 9.856169820212966E-01 MA= 3.572170843984495E+02 TA= 3.571221738148128E+02
 A = 9.999947139070439E-01 AD= 1.016698244889530E+00 PR= 3.652534468934519E+02''',
									  needSystemGM=False)
	return jdtomjd(jd[0]),E[0]

def EMB_xyz_at_times(times):
	t0,E = EMB_ephem()
	(a,e,I,Omega,pomega,M0,nil) = E
	GM = 2.9591310798672560E-04 #AU^3/d^2
	dMdt = 	sqrt(GM / a**3)
	obsxyz = []
	for t in times:
		M = M0 + dMdt * (t - t0)
		(x,v) = phase_space_coordinates_from_orbital_elements(a,e,I,Omega,pomega,M,GM)
		obsxyz.append(x)
	obsxyz = array(obsxyz)
	return obsxyz

def plot_ellipses(data):
	allra  = data[0]
	alldec = data[1]
	allrad = data[2]
	I = argsort(-allrad)
	for ra,dec,rad in zip(allra[I],alldec[I],allrad[I]):
		facealpha = min(1.0, max(0.005, 0.1/(rad**2)))
		edgealpha = min(1.0, max(0.05, 0.1/rad))
		gca().add_artist(Ellipse([ra,dec], 2.*rad/cosdeg(dec), 2.*rad,
					 fc='k', ec='none', alpha=facealpha))
		# under our shitty Py 2.6 build edges aren't ever transparent,
		# so apply edgealpha as a width not an alpha
		gca().add_artist(Ellipse([ra,dec], 2.*rad/cosdeg(dec), 2.*rad,
					 fc='none', ec='k', lw=edgealpha))
		# Ideally:
		#gca().add_artist(Ellipse([ra,dec], 2.*rad/cosdeg(dec), 2.*rad,
		#			 fc='none', ec='k', alpha=edgealpha, lw=0.3))

class CometMCMCxv:
	#             AU,  AU/day
	paramnames = ['x', 'v']
	def __init__(self):
		self.x = zeros(3)
		self.v = zeros(3)
		self.xstep = 0.0003
		self.vstep = 0.00004
		self.reset_counters()
		# epoch of the orbital elements (ie, M)
		self.epoch = None
		# Step size for comet trajectory
		self.dtdays = 2.
		self.times = None
		self.earthxyz = None

	def set_true_params(self):
		# from test-holmes-1.txt
		s = '''2454418.500000000 = A.D. 2007-Nov-14 00:00:00.0000 (CT)
   1.248613009072901E+00  2.025389080777020E+00  8.242272661173784E-01
  -7.973747689948190E-03  9.391385050634651E-03  1.214741189900433E-03
   1.454305558379576E-02  2.518052017168866E+00  3.997656275386785E-03
   '''
		x,v,jd = jpl.parse_phase_space(s)
		x = x[0]
		v = v[0]
		jd = jd[0]
		self.set_x(x)
		self.set_v(v)
		self.set_epoch(jdtomjd(jd))

	def get_params(self):
		return hstack((self.get_x(), self.get_v()))

	def get_x(self):
		return self.x.copy()

	def get_v(self):
		return self.v.copy()

	def set_x(self, x):
		self.x[:] = x[:]

	def set_v(self, v):
		self.v[:] = v[:]

	def set_params(self, p):
		self.x[:] = p[:3]
		self.v[:] = p[3:]

	def propose_params(self):
		i = numpy.random.randint(2)
		x0 = self.get_x()
		v0 = self.get_v()
		if i == 0:
			x0 += self.xstep * numpy.random.normal(size=3)
		else:
			v0 += self.vstep * numpy.random.normal(size=3)
		self.proposed_param = i
		return hstack((x0, v0))

	def reset_counters(self):
		self.nproposed = zeros(2)
		self.naccepted = zeros(2)

	def tally(self, accepted, linknumber):
		self.nproposed[self.proposed_param] += 1
		if accepted:
			self.naccepted[self.proposed_param] += 1
		if linknumber % 10 == 0:
			print 'link:', linknumber, 'lnp', self.lnp,
			print ' -- ratios:', '[' + ', '.join(['%.3f' % x for x in (self.naccepted.astype(float)/self.nproposed)]) + ']'
			print 'params:',
			#for p,n in zip(self.get_params(), CometMCMCxv.paramnames):
			#	print n, '%.5f' % p,
			x = self.get_x()
			v = self.get_v()
			print 'x [ %.5f, %.5f, %.5f ]' % (x[0],x[1],x[2]),
			print 'v [ %.5f, %.5f, %.5f ]' % (v[0],v[1],v[2])
			print

		#if (linknumber+1) % self.plotinterval == 0:
		#	print 'Sample', self.samplenumber, 'link', linknumber, 'plotting'
		#	'''
		#	figure(1)
		#	self.plot_radecs()
		#	#figure(2)
		#	#self.plot_xyzs()
		#
		#	if self.saveplots:
		#		print 'Sample', self.samplenumber, 'link', linknumber, 'saving'
		#		figure(1)
		#		savefig('mcmc-%05i.png' % linknumber)
		#		#figure(2)
		#		#savefig('xyz-%05i.png' % linknumber)
		#		'''			
	def plot_radecs(self):
		ras,decs = self.radec_at_times()
		if any(isinf(ras)) or any(isinf(decs)):
			print 'Infinite RA,Dec'
		if not(all(isfinite(ras)) and all(isfinite(decs))):
			print 'non-finite RA,Dec'
		print 'RA range', min(ras), max(ras)
		print 'Dec range', min(decs), max(decs)
		plot(ras, decs, 'r-', alpha=0.25)
		axis([20,80,30,60])

	'''
	def plot_xyzs(self):
		(xhat,yhat,zhat) = orbital_vectors_from_orbital_elements(self.get_i_rad(), self.get_Om_rad(), self.get_om_rad())
		M = linspace(0, 2*pi, 360)
		C = array([position_from_orbital_vectors(xhat,yhat, self.get_a(), self.get_e(), mi)
				   for mi in M])
		subplot(2,2,1)
		plot(C[:,0], C[:,1], 'g-', alpha=0.1)
		axis('equal')
		subplot(2,2,3)
		plot(C[:,0], C[:,2], 'g-', alpha=0.1)
		axis('equal')
		subplot(2,2,2)
		plot(C[:,2], C[:,1], 'g-', alpha=0.1)
		axis('equal')
		'''

	def lnposterior(self, data):
		lnp = self.lnprior() + self.lnlikelihood(data)
		self.lnp = lnp
		return lnp

	def lnprior(self):
		return 0

	def set_epoch(self, mjd):
		self.epoch = mjd

	def set_date_range(self, tmin, tmax):
		self.times = arange(tmin, tmax, self.dtdays)
		self.earthxyz = EMB_xyz_at_times(self.times)

	def get_orbital_elements(self):
		return orbital_elements_from_phase_space_coordinates(self.get_x(), self.get_v(), GM_sun)

	def radec_at_times(self, times=None, light_travel=True):
		cometras = []
		cometdecs = []
		C = list(self.get_orbital_elements()) + [GM_sun]
		CM0 = C[5]
		CdMdt = sqrt(GM_sun / C[0]**3)

		if times is None:
			times = self.times
		if times is self.times:
			earthxyz = self.earthxyz
		else:
			earthxyz = EMB_xyz_at_times(times)

		for ex,t in zip(earthxyz, times):
			C[5] = CM0 + CdMdt * (t - self.epoch)
			(r,d) = orbital_elements_to_radec(C, ex)
			cometras.append(r)
			cometdecs.append(d)
		return (array(cometras), array(cometdecs))

	def lnlikelihood(self, data):
		(ras,decs,radii) = data
		imgxyz = radectoxyz(ras, decs)
		imgrad = sqrt(deg2distsq(radii)) * 0.5

		gQ = 10. * 4.*pi/(imgrad**2) / (max(self.times)-min(self.times))

		## HACK -- just convert to orbital elements and carry on as usual.
		E = self.get_orbital_elements()
		a = E[0]
		e = E[1]
		I_rad = E[2]
		Omega_rad = E[3]
		pomega_rad = E[4]
		M_rad = E[5]

		lnl = comet_lnlikelihood_gaussian(a, e, I_rad, Omega_rad, pomega_rad, M_rad,
						  self.epoch, self.times, self.earthxyz,
						  imgxyz, imgrad, gQ)
		return lnl


def best_and_secondbest(X):
	ibest = argmin(X)
	if ibest == 0:
		i2 = argmin(X[ibest+1:]) + ibest+1
		return (ibest, i2)
	if ibest == (len(X)-1):
		i1 = argmin(X[:ibest])
		return (ibest, i1)
	i1 = argmin(X[:ibest])
	i2 = argmin(X[ibest+1:]) + ibest+1
	if X[i1] < X[i2]:
		isecond = i1
	else:
		isecond = i2
	return (ibest, isecond)

# Read the given data pickle, initialize a number of chains and return them.
def init_chains(imgdata, nchains):
	data = (array(imgdata['center-ras']),
			array(imgdata['center-decs']),
			array(imgdata['radii']))

	dates = imgdata['dates']
	print 'Number of images:', len(dates)

	mjds = array([datetomjd(d) if d is not None else 0 for d in dates])
	I = (mjds > 0)

	tmin = min(mjds[I]) - 30
	tmax = max(mjds[I]) + 30

	print 'Date range: MJD', tmin, tmax
	print 'Date range:', mjdtodate(tmin), mjdtodate(tmax)

	# Initialize orbit using images close to the median date.
	medmjd = median(mjds[I])
	print '%i valid dates' % sum(I), 'median', medmjd
	nearby = (abs(mjds - medmjd) < 7)
	print 'number of images within 7 days:', sum(nearby)

	epoch = medmjd
	# Move the epoch to the nearest integer MJD.
	epoch = float(round(epoch))

	nra  = data[0][nearby]
	ndec = data[1][nearby]
	nrad = data[2][nearby]
	nt = mjds[nearby]
	N = len(nra)
	print 'N', N
	#print 'nra,ndec,nrad,nt', nra.shape, ndec.shape, nrad.shape, nt.shape
	sqrtW = 1./nrad
	#print 'sqrtW', sqrtW.shape

	X = zeros((N,2))
	X[:,0] = 1. * sqrtW
	X[:,1] = (nt-epoch) * sqrtW
	y = zeros((N,2))
	y[:,0] = nra * sqrtW
	y[:,1] = ndec * sqrtW
	(b,resid,rank,eigs) = lstsq(X, y)
	print 'b', b
	ra0,dec0 = b[0,:]
	dradt,ddecdt = b[1,:]

	clf()
	subplot(2,1,1)
	errorbar(nt-epoch, nra, yerr=nrad, fmt=None)
	a = axis()
	tt = array(a[:2])
	print 'tt', tt
	plot(tt, ra0 + dradt * tt, 'b-')
	xlim(a[0],a[1])
	ylim(ra0-5, ra0+5)
	ylabel('RA (deg)')
	subplot(2,1,2)
	errorbar(nt-epoch, ndec, yerr=nrad, fmt=None)
	plot(tt, dec0 + ddecdt * tt, 'b-')
	xlim(a[0],a[1])
	ylim(dec0-5, dec0+5)
	ylabel('Dec (deg)')
	xlabel('dt (days)')
	savefig('ravst.png')

	allra  = data[0]
	alldec = data[1]
	allrad = data[2]

	figure(1)
	clf()
	for ra,dec,rad in zip(allra,alldec,allrad):
		gca().add_artist(Ellipse([ra,dec], 2.*rad/cosdeg(dec), 2.*rad,
								 fc='0.5', ec='0.5', alpha=0.01))
	for ra,dec,rad in zip(nra,ndec,nrad):
		gca().add_artist(Ellipse([ra,dec], 2.*rad/cosdeg(dec), 2.*rad,
								 fc='b', ec='b', alpha=0.01))
	tt = arange(-7, 8)
	plot(ra0 + tt * dradt, dec0 + tt * ddecdt, 'r-')
	axis([30,70,30,70])
	savefig('nearby.png')

	# Find velocity and compute orbital elements from x,v.
	dt = 1.
	ra1  = ra0  + dradt  * dt
	dec1 = dec0 + ddecdt * dt
	xyz0 = radectoxyz(ra0, dec0)[0]
	xyz1 = radectoxyz(ra1, dec1)[0]
	print 'xyz:', xyz0, xyz1
	# xyz are in celestial coords; convert to solar-system coords
	(antieq, antisol, antipole) = ecliptic_basis(eclipticangle = -23.439281)
	xyz0 = xyz0[0] * antieq + xyz0[1] * antisol + xyz0[2] * antipole
	xyz1 = xyz1[0] * antieq + xyz1[1] * antisol + xyz1[2] * antipole

	# Find observer (EMB) positions at epoch
	t0_emb,E_emb = EMB_ephem()
	E_emb = list(E_emb)
	(a,e,I,Omega,pomega,Mx,nil) = E_emb
	dMdt_emb = sqrt(GM_sun / a**3)
	M0 = Mx + dMdt_emb * (epoch - t0_emb)
	(embx0,v) = phase_space_coordinates_from_orbital_elements(a,e,I,Omega,pomega,M0,GM_sun)
	M1 = Mx + dMdt_emb * ((epoch+dt) - t0_emb)
	(embx1,v) = phase_space_coordinates_from_orbital_elements(a,e,I,Omega,pomega,M1,GM_sun)
	#print 'xyz0,xyz1', xyz0, xyz1
	#print 'embx0,1', embx0, embx1

	ctrue = CometMCMCxv()
	ctrue.set_true_params()
	ctrue.set_date_range(tmin, tmax)

	figure(1)
	clf()
	plot_ellipses(data)
	ras,decs = ctrue.radec_at_times()
	plot(ras, decs, 'b-')

	rads = data[2].copy()
	rads.sort()
	print 'Smallest 10 images:', rads[:10], 'deg'
	print 'Largest angular motion in %g days:' % ctrue.dtdays
	btw = [degrees_between(ras[i-1], decs[i-1], ras[i], decs[i])
		   for i in range(1, len(ras))]
	print max(btw), 'deg'

	'''
	figure(2)
	clf()
	(xhat,yhat,zhat) = orbital_vectors_from_orbital_elements(ctrue.get_i_rad(), ctrue.get_Om_rad(), ctrue.get_om_rad())
	M = linspace(0, 2*pi, 360)
	xyzs = array([position_from_orbital_vectors(xhat,yhat, ctrue.get_a(), ctrue.get_e(), mi)
				  for mi in M])
	print 'xyzs shape:', xyzs.shape
	subplot(2,2,1)
	E = ctrue.earthxyz
	C = xyzs
	plot(C[:,0], C[:,1], 'r-')
	plot(E[:,0], E[:,1], 'b-')
	plot([0],[0], 'yo')
	axis('equal')
	xlabel('x')
	ylabel('y')
	subplot(2,2,3)
	plot(C[:,0], C[:,2], 'r-')
	plot(E[:,0], E[:,2], 'b-')
	plot([0],[0], 'yo')
	axis('equal')
	xlabel('x')
	ylabel('z')
	subplot(2,2,2)
	plot(C[:,2], C[:,1], 'r-')
	plot(E[:,2], E[:,1], 'b-')
	plot([0],[0], 'yo')
	axis('equal')
	xlabel('z')
	ylabel('y')
	savefig('truexyz.png')
	'''

	print 'True params:', ctrue.get_params()
	print 'at epoch', ctrue.epoch
	print 'lnlikelihood of true parameters:', ctrue.lnlikelihood(data)
	print

	mcmcs = []
	for nc in range(nchains):
		while True:
			#R = 10.0**numpy.random.uniform(-1.5,-0.5)
			R = 10.0**(linspace(-0.5, 0.5, nchains)[nc])
			# Comet position and velocity...
			cx0 = embx0 + R * xyz0
			cx1 = embx1 + R * xyz1
			cv0 = (cx1 - cx0) / dt
			print 'R = ', R
			print 'angle between cx0, cv0:', rad2deg(arccos(dot(cx0, cv0) / norm(cx0) / norm(cv0)))
			print 'kinetic energy:', 0.5 * dot(cv0,cv0)
			print 'potential energy:', GM_sun / norm(cx0)
			E = 0.5 * dot(cv0,cv0) - GM_sun / norm(cx0)
			# FIXME -- with fixed radii, this could loop infinitely...
			if E < 0:
				print 'Keeping this one.'
				break

		# Start MCMC
		c = CometMCMCxv()
		c.set_epoch(epoch)
		c.set_date_range(tmin, tmax)
		c.set_x(cx0)
		c.set_v(cv0)
		c.R = R
		print 'Initial params:', c.get_params()
		print 'lnlikelihood of initial parameters:', c.lnlikelihood(data)
		print
		print 'Chose epoch:', epoch

		c.chain = []
		c.bestlnp = -1e100
		c.bestparams = None
		c.samplenumber = nc + 1

		figure(1)
		c.plot_radecs()

		mcmcs.append(c)

	print 'Chose radii:', [c.R for c in mcmcs]

	#figure(2)
	#c.plot_xyzs()

	figure(1)
	savefig('mcmc-00000.png')

	return {'chains': mcmcs,
		'ctrue': ctrue,
		'data': data,
		'dates': dates,
		'tmin':tmin,
		'tmax':tmax,
		'epoch':epoch}


def metropolis_thread_target(data, c, ns, step, beta, rq=None):
	#print 'Args:', data, c, ns, step, beta, rq
	(thislnp, thisparams, thischain) = metropolis(data, c, ns, startlink=step, beta=beta)
	c.chain += thischain
	if thislnp > c.bestlnp:
		c.bestlnp = thislnp
		c.bestparams = thisparams
	if rq:
		rq.put(c)
	print 'worker best lnp:', c.bestlnp
	return c

def metropolis_proxy(X):
	return metropolis_thread_target(*X)

def run_chain(D, nsteps, beta, prefix):
	mcmcs = D['chains']
	ctrue = D['ctrue']
	data = D['data']
	dates = D['dates']
	tmin = D['tmin']
	tmax = D['tmax']
	epoch = D['epoch']

	figure(1)
	clf()
	plot_ellipses(data)
	ras,decs = ctrue.radec_at_times()
	plot(ras, decs, 'b-')

	ctemp = CometMCMCxv()

	if multiprocessing:
		workers = Pool(processes = 8)
	else:
		workers = None

	step = 0
	snapshot = 500
	while step < nsteps:
		ns = min(nsteps, step + snapshot)

		if multiprocessing:
			print 'Distributing jobs to workers...'
			args = []
			for c in mcmcs:
				args.append((data, c, ns, step, beta))
			pre = [c.bestlnp for c in mcmcs]
			mcmcs = workers.map(metropolis_proxy, args, 1)
			post = [c.bestlnp for c in mcmcs]
			print 'Workers are finished!'
			print 'pre best lnp:', pre
			print 'post best lnp:', post

		else:
			threads = []
			# With 'threading', we don't really need the queues!
			qs = []
			print 'Launching threads...'
			for ic,c in enumerate(mcmcs):
				# Per-process result queue.
				q = Queue()
				t = Thread(target=metropolis_thread_target, args=(data, c, ns, step, beta, q), name='chain-%i' % ic)
				t.start()
				threads.append(t)
				qs.append(q)
			print 'Gathering threads...'
			for ic,t in enumerate(threads):
				mcmcs[ic] = qs[ic].get()
				t.join()
			print 'Gathered threads!'

		for c in mcmcs:
			# Pull out best params, convert to orbital elements
			if c.bestparams is None:
				print 'comet', c, 'has no best params!'
				continue
			ctemp.set_epoch(c.epoch)
			ctemp.set_params(c.bestparams)
			ctemp.times = c.times
			ctemp.earthxyz = c.earthxyz
			print 'Best params thus far:', c.bestparams
			print 'Best lnprob:', c.bestlnp
			print 'Orbital elements:', ctemp.get_orbital_elements()

			figure(1)
			c.plot_radecs()

		figure(1)
		savefig('%s-%05i.png' % (prefix, ns))

		step = ns

	return {'chains': mcmcs,
		'ctrue': ctrue,
		'data': data,
		'dates': dates,
		'tmin':tmin,
		'tmax':tmax,
		'epoch':epoch}

	if False:
		# Remove the radial component of velocity...
		vr = xyz0 * dot(cv0, xyz0) / norm(xyz0)**2
		cv0 -= vr

def plot_from_chain(chainpicklefn):
	p = unpickle_from_file(chainpicklefn)
	chain = p['chain']
	tmin = p['tmin']
	tmax = p['tmax']
	data = p['data']
	epoch = p['epoch']
	(data_ras, data_decs, data_radii) = data

	lnps = array([lnp for (lnp,p) in chain])
	params = [p for (lnp,p) in chain]

	if True:
		clf()
		plot(lnps, 'k.')
		xlabel('iteration')
		ylabel('lnp')
		savefig('iter-lnp.png')

		clf()
		hist(exp(lnps - median(lnps)), 50)
		xlabel('p')
		savefig('hist-p.png')

		for i in range(len(params[0])):
			parami = [p[i] for p in params]
			clf()
			plot(parami, 'k.')
			xlabel('iteration')
			ylabel(CometMCMC.paramnames[i])
			savefig('iter-%s.png' % CometMCMC.paramnames[i])

			clf()
			hist(parami, 50)
			xlabel(CometMCMC.paramnames[i])
			savefig('hist-%s.png' % CometMCMC.paramnames[i])

	# Trajectory:
	ibest = argmax(lnps)
	pbest = params[ibest]
	lnpbest = lnps[ibest]

	c = CometMCMC()
	c.set_date_range(tmin, tmax)

	print 'getting orbit for true params...'
	c.set_true_params()
	(true_ras, true_decs) = c.radec_at_times()

	c.set_epoch(epoch)
	c.set_date_range(tmin, tmax)
	clf()
	for i in range(1):
		# sample from second half of chain.
		r = numpy.random.randint(len(params)/2) + len(params)/2
		print 'plotting orbit for params:', params[r]
		c.set_params(params[r])
		(ras,decs) = c.radec_at_times()
		p1 = plot(ras, decs, 'k-', alpha=0.25, linewidth=2)

	p2 = plot(true_ras, true_decs, 'r.-', alpha=1)

	# label some dates
	ltimes = [datetime.datetime(x/12, x%12+1, 1) for x in range(2007*12+7, 2008*12+3)]
	(lras,ldecs) = c.radec_at_times(array([datetomjd(d) for d in ltimes]))

	c.set_params(pbest)
	lnp = c.lnposterior(data)
	assert(lnp == lnpbest)
	(ras,decs) = c.radec_at_times()
	p3 = plot(ras, decs, 'b-', alpha=1)

	print '%i images' % (len(data_ras))
	# a = gca()
	for ra,dec,radius in zip(data_ras,data_decs,data_radii):
		gca().add_artist(Ellipse([ra,dec], 2.*radius/cosdeg(dec), 2.*radius,
								 fc='none', ec='k',
								 #fc='k', ec='none',
								 alpha=0.1))
	plot(lras, ldecs, 'ko')
	for t,ra,dec in zip(ltimes,lras,ldecs):
		text(ra+0.5, dec+0.5, str(t.date()), horizontalalignment='left', verticalalignment='bottom',
		     bbox=dict(facecolor='w', alpha=0.25))

	legend((p2,p1,p3), ("``Truth''", 'Samples', 'Best sample'))
	xlabel('RA (deg)')
	ylabel('Dec (deg)')
	axis([30,70,30,60])
	savefig('trajectory.png')

	axis([0,180,0,60])
	savefig('trajectory2.png')

if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) < 1:
		print 'Usage: %s <command>' % sys.argv[0]
		print 'Commands include:'
		print '   init'
		print '   burn'
		print '   run'
		print '   resample'
		sys.exit(-1)
	cmd = args[0]
	cmdargs = args[1:]
	if cmd == 'init':
		Nchains = 0
		if len(cmdargs):
			Nchains = int(cmdargs[0])
		if not Nchains:
			Nchains = 8
		print 'Initializing', Nchains, 'chains'
		imgdata = unpickle_from_file('imagewcs.pickle')
		D = init_chains(imgdata, Nchains)
		pickle_to_file(D, 'init.pickle')
	elif cmd == 'burn':
		D = unpickle_from_file('init.pickle')
		for beta in arange(0.1,1.01,0.2):
			D = run_chain(D, 200, beta, 'burn-%04.2f' % beta)
		pickle_to_file(D, 'run.pickle')
	elif cmd == 'run':
		D = unpickle_from_file('run.pickle')
		D = run_chain(D, 1000, 1.0, 'mcmc')
		pickle_to_file(D, 'run.pickle')
	elif cmd == 'resample':
		Nchains = 0
		if len(cmdargs):
			Nchains = int(cmdargs[0])
		if not Nchains:
			Nchains = 8
		print 'Resampling', Nchains, 'chains'

		D = unpickle_from_file('run.pickle')
		chains = D['chains']
		tmin = D['tmin']
		tmax = D['tmax']

		maxLs = []
		rLs = []
		for c in chains:
			L = []
			for (lnp,params) in c.chain:
				L.append(lnp)
			L = array(L)
			maxL = max(L)
			maxLs.append(maxL)
			print 'max L:', maxL
			rLs.append(sum(exp(L - maxL)))
			print 'rL:', sum(exp(L-maxL))
		maxLs = array(maxLs)
		maxLs -= max(maxLs)
		T = exp(maxLs) * array(rLs)
		print 'Total probabilities:', T
		T /= sum(T)
		cT = cumsum(T)

		samples = []
		for i in range(Nchains):
			u = numpy.random.uniform()
			# which chain do we pull it from?
			I = argmax(cT > u)
			c = chains[I]
			# sample a link from this chain.
			I = numpy.random.randint(len(c.chain))
			(lnp,params) = c.chain[I]

			newc = CometMCMCxv()
			newc.epoch = c.epoch
			newc.R = c.R
			newc.set_date_range(tmin, tmax)
			newc.chain = []
			newc.bestlnp = -1e100
			newc.bestparams = None
			newc.samplenumber = i+1
			newc.set_params(params)

			samples.append(newc)

		D['chains'] = samples
		pickle_to_file(D, 'run.pickle')

	elif cmd == 'plots':
		D = unpickle_from_file('run.pickle')

		Rs = []
		Ls = []
		for c in D['chains']:
			R = []
			L = []
			for (lnp,params) in c.chain:
				R.append(norm(params[:3]))
				L.append(lnp)
			Rs.append(R)
			Ls.append(L)
		clf()
		for R in Rs:
			plot(R)
		xlabel('MCMC step')
		ylabel('Comet orbital radius at epoch')
		savefig('R.png')

		clf()
		for L in Ls:
			plot(L)
		xlabel('MCMC step')
		ylabel('log-prob of orbital parameters')
		savefig('lnp.png')
		
	else:
		print 'Unknown command', cmd
	sys.exit(0)

	#import cProfile
	cfn = 'chain.pickle'
	#cProfile.run("main('imagewcs.pickle', cfn, 100)")
	#main('imagewcs.pickle', cfn, 100)
	#sys.exit(0)
	
	if not os.path.exists(cfn):
		main('imagewcs.pickle', cfn, 1000000)
	#print 'bailing out'
	#sys.exit(0)
	#print 'Reading from', cfn
	#plot_from_chain(cfn)

