import matplotlib
if __name__ == '__main__':
	matplotlib.use('Agg')
import os
import datetime
import time
import sys
from math import pi,sqrt
import numpy as np
import pylab as plt
from matplotlib.patches import Polygon
from astrometry.util import EXIF
from astrometry.util.starutil_numpy import jdtomjd, radectoxyz, deg2distsq, datetomjd, mjdtodate, arcsec2distsq, xyztoradec, arcsec_between
from astrometry.util.file import *
from astrometry.util.util import Tan
from mcmc_comet import CometMCMCxv

pfn = 'exifplots.pickle'

def tally(S):
	u,I = np.unique(np.array(S), return_inverse=True)
	#print 'Makes:', u
	N = []
	for i in range(len(u)):
		n = np.sum(i == I)
		N.append(n)
		#print u[i], n
	I = np.argsort(-np.array(N))
	return [(u[i], N[i]) for i in I]

def exifplots(C, data, exifs, makeplots=None):
	if makeplots is None:
		makeplots = range(1,10+1)

	print 'Pickling', pfn
	pickle_to_file((C,data,exifs), pfn)

	cams = []
	makes = []
	softs = []
	models = []

	allkeys = set()
	for jpegfn,exif in exifs:
		if not len(exif.keys()):
			continue

		make = str(exif.get('Image Make', ''))
		mod  = str(exif.get('Image Model', ''))
		cams.append(make + ' / ' + mod)
		mw = make.split()
		man = len(mw)>0 and mw[0] or ''
		mod = ' '.join(mw[1:]) + ' ' + mod
		models.append((man, mod))

		makes.append(man)#str(exif.get('Image Make', '')))
		softs.append(str(exif.get('Image Software', '')))

		print
		print jpegfn
		allkeys.update(exif.keys())
		for k in ['Image Software',
				  'Image Model',
				  'Image Make',
				  'GPS GPSVersionID',
				  'GPSDateStamp',
				  'GPSTimeStamp',
				  'Image Artist',
				  'Image BitsPerSample',
				  'Image Copyright',
				  'Image GPSInfo',
				  'Image ImageDescription',
				  'Image PhotometricInterpretation',
				  'EXIF ShutterSpeedValue',
				  'EXIF MeteringMode',
				  'EXIF ApertureValue',
				  'EXIF FNumber',
				  'EXIF FocalLength',
				  'EXIF ExposureTime',
				  'EXIF ISOSpeedRatings',
				  'EXIF FocalLengthIn35mmFilm',
				  'EXIF FocalPlaneResolutionUnit',
				  'EXIF FocalPlaneXResolution',
				  'EXIF FocalPlaneYResolution']:
			if k in exif:
				print k, exif[k]

	allkeys = list(allkeys)
	allkeys.sort()
	print 'All EXIF keys seen:', allkeys

	print
	print 'Makes + Models:'
	cams = [x for x in cams if len(x) > 1]
	for u,n in tally(cams):
		print '  ', n, ' ', u

	makes = [x for x in makes if len(x)]
	print 'Manufacturers:'
	for u,n in tally([m.split()[0].upper() for m in makes]):
		print '  ', n, ' ', u

	softs = [x for x in softs if len(x)]
	print 'Software:'
	for u,n in tally(softs):
		print '  ', n, ' ', u

	matplotlib.rc('font', family='computer modern roman')
	matplotlib.rc('font', size=14)
	matplotlib.rc('text', usetex=True)
	plt.figure(figsize=(6,6))
	spargs = dict(bottom=0.15, left=0.15, right=0.95, top=0.95)

	formats = ['.pdf', '.png', '.eps']

	prefix = 'exif1'
	if 1 in makeplots:
		plt.clf()
		plt.subplots_adjust(**spargs)
		#mn,mx = 54290, 54550
		# mn,mx = C.tmin, C.tmax
		mn,mx = datetomjd(datetime.datetime(2007,7,1)),datetomjd(datetime.datetime(2008,4,1))
		D = np.array([d for (wcs,d) in data if d > 0])
		print 'Total # data:', len(data)
		print 'Number with mjd > 0:', len(D)
		(n, b, p) = plt.hist(D, range=(mn,mx), bins=50, histtype='step', color='k')
		print 'Number within histogram range:', sum(n)
		ltimes = [datetime.datetime(x/12, x%12+1, 1)
				  for x in range(2007*12 + 5, 2008*12 + 4, 1)]
		plt.xticks([datetomjd(d) for d in ltimes],
				   [str(d.date()) if ((i%2) == 0) else '' for i,d in enumerate(ltimes)],
				   rotation=15, verticalalignment='top', horizontalalignment='right')
		plt.xlabel('date')
		plt.axvline(datetomjd(datetime.datetime(2007, 10, 24)),
					color='r', linestyle='--', lw=2)
		plt.ylabel('number of images')
		plt.xlim(mn,mx)
		plt.axhline(0., color='0.75')
		plt.ylim(-0.1 * np.max(n), 1.1 * np.max(n))
		for f in formats:
			plt.savefig(prefix + f)

	# manufacturers
	cleanmakes = [m.split()[0].upper() for m in makes]
	X = tally(cleanmakes)
	man = np.array([x[0] for x in X])
	nman = np.array([x[1] for x in X])
	ncut = 4
	I = (nman < 4)
	J = (nman >= 4)
	man  = np.append(man[J], ['Other'])
	nman = np.append(nman[J], [sum(nman[I])])

	print 'man', man
	print 'nman', nman

	print 'Total # of images with Manufacturer EXIF tags:', sum(nman)

	prefix = 'exif2'
	plotfn = prefix+'.pdf'
	#if not os.path.exists(plotfn):
	if 2 in makeplots:
		W = 0.8
		plt.clf()
		plt.subplots_adjust(**spargs)
		X = np.arange(len(man))
		plt.bar(X, nman, color='none', ec='k', width=W)
		plt.xticks([])
		for x,n,m in zip(X+0.5, nman, man):
			plt.text(x-0.2, n, m, ha='left', va='bottom', rotation=45)

		for x,n,m in zip(X+0.5, nman, man):
			print 'Manufacturer', m

			# 300D = Rebel
			# 350D = Rebel XT
			# 400D = XTi
			# 450D = XSi
			modmap = {
				'REBEL': '300D',
				'REBEL XT': '350D',
				'REBEL XTi': '400D',
				}
			mods = [mod for mk,mod in models if len(mk)>0 and mk.split()[0].upper() == m]
			# print 'mods:'
			# for u,ni in tally(mods):
			# 	print '  ', ni, ' ', u
			mods = [mod.replace('Canon ','').replace('DIGITAL ','').replace('EOS ','').replace('DIGITAL','') for mod in mods]
			mods = [mod.replace('NIKON ','').replace('CORPORATION ','') for mod in mods]
			# olympus
			mods = [mod.replace('OPTICAL CO.,LTD ','') for mod in mods]
			# print 'mods:'
			# for u,ni in tally(mods):
			# 	print '  ', ni, ' ', u
			mods = [mod.strip() for mod in mods if len(mod.strip()) > 0]
			mods = [modmap.get(mod,mod) for mod in mods]
			X = tally(mods)
			print 'Models:'
			for u,ni in X:
				print '  ', ni, ' ', u

			top = n
			for u,ni in X:
				if ni < 5:
					break
				if len(u) > 4:
					u = u[:4]+'..'
				plt.text(x-0.5+W/2., top-(ni/2.), u, ha='center', va='center', rotation=0, fontsize=8, color='0.5') #alpha=0.5)
				plt.plot([x-0.5, x-0.5+W], [top-ni,top-ni], '-', color='0.5')#color='k', alpha=0.5)
				top -= ni
			



		plt.xlabel('manufacturer')
		plt.ylabel('number of images')
		margin = 1.2
		plt.xlim(W-margin,len(man)-1+margin)
		plt.axhline(0., color='k', alpha=0.25)
		# note magic to make everything line up
		plt.ylim(-0.1 * (1.2 / 1.1) * np.max(nman), 1.2 * np.max(nman))
		for f in formats:
			plt.savefig(prefix + f)


	# Image resolution and extent
	scale = []
	area = []
	for wcs,d in data:
		cd = wcs[4:8]
		pixscale = sqrt(abs(cd[0]*cd[3] - cd[1]*cd[2]))
		W,H = wcs[8],wcs[9]
		scale.append(pixscale)
		area.append(pixscale**2 * W * H)
	area = np.array(area)
	pixscale = np.array(scale)

	prefix = 'exif3'
	plotfn = prefix+'.pdf'
	if 3 in makeplots:
		plt.clf()
		(n, b, p) = plt.hist(np.log10(3600*pixscale), 40, histtype='step', color='k')
		plt.xlabel('image pixel scale (arcsec/pixel)')
		plt.xticks([-1, 0, 1, 2, 3], ['0.1', '1', '10', '100', '1000'])
		plt.ylabel('Number of images')
		plt.axhline(0., color='k', alpha=0.25)
		plt.ylim(-0.1 * np.max(n), 1.1 * np.max(n))
		for f in formats:
			plt.savefig(prefix + f)


	prefix = 'exif4'
	if 4 in makeplots:
		plt.clf()
		plt.subplots_adjust(**spargs)
		(n, b, p) = plt.hist(np.log10(area), 40, histtype='step', color='k')
		plt.axhline(0., color='0.75') #color='k', alpha=0.25) #color='0.25')
		plt.xlabel('image area (deg$^2$)')
		plt.xticks([-2, -1, 0, 1, 2, 3],
				   [0.01, 0.1, 1, 10, 100, 1000],
				   fontsize='medium')
		plt.ylabel('number of images')
		plt.ylim(-0.1 * np.max(n), 1.1 * np.max(n))
		for f in formats:
			plt.savefig(prefix + f)

	prefix = 'exif5'
	if 5 in makeplots:
		plt.clf()
		(n, b, p) = plt.hist(np.log10(np.sqrt(area)), 40, histtype='step', color='k')
		plt.xlabel('image extent (deg)')
		plt.ylabel('number of images')
		plt.axhline(0., color='k', alpha=0.25)
		plt.ylim(-0.1 * np.max(n), 1.1 * np.max(n))
		for f in formats:
			plt.savefig(prefix + f)

	def Rtofloat(self):
		#print self.num, '/', self.den
		return float(self.num) / float(self.den)
	EXIF.Ratio.__float__ = Rtofloat

	prefix = 'exif6'
	if 6 in makeplots:
		K = 'EXIF ExposureTime'
		exptime = []
		for fn,exif in exifs:
			if K in exif:
				exptime.append(float(exif.get(K).values[0]))
		exptime = np.array(exptime)
		exptime = exptime[exptime > 0]
		print len(exptime), 'images have', K
		plt.clf()
		plt.subplots_adjust(**spargs)
		(n, b, p) = plt.hist(np.log10(exptime), 20, histtype='step', color='k')
		plt.xlabel('image exposure time (s)')
		plt.xticks([-1, 0, 1, 2, 3],
				   [0.1, 1, 10, 100, 1000])
		plt.ylabel('number of images')
		plt.axhline(0., color='0.75') #color='k', alpha=0.25)
		plt.ylim(-0.1 * np.max(n), 1.1 * np.max(n))
		for f in formats:
			plt.savefig(prefix + f)

	pexif = 0.71
	DD = datetime.datetime(2007, 10, 24, 23, 11, 13)
	D = datetomjd(DD)
	H = C.find_times_slice_within_12_hours_of_mjd(D)
	assert(H.sum() > 0)
	ptime = (1. - pexif) * C.ptime
	# MAGIC: the following lines rely on dtdays being in days!
	ptime[H] += pexif
	ptime /= (np.sum(ptime) * C.dtdays)

	print 'sum ptimes:', sum(C.ptime), sum(ptime)
	mn,mx = C.tmin, C.tmax

	prefix = 'exif7'
	if 7 in makeplots:
		plt.clf()
		plt.subplots_adjust(**spargs)
		p1 = plt.plot(C.times, C.ptime, 'k-', alpha=0.5, lw=2)
		p2 = plt.plot(C.times, ptime, 'k-')

		ltimes = [datetime.datetime(x/12, x%12+1, 1)
				  for x in range(2007*12 + 5, 2008*12 + 4, 1)]
		plt.xticks([datetomjd(d) for d in ltimes],
				   [str(d.date()) if ((i%2) == 0) else ''
					for i,d in enumerate(ltimes)],
				   rotation=15, verticalalignment='top',
				   horizontalalignment='right')
		plt.xlabel('date')
		#plt.axvline(datetomjd(datetime.datetime(2007, 10, 24)),
		#			color='r', linestyle='--', lw=2)
		plt.ylabel('probability density (d$^{-1}$)')
		plt.yticks([])
		plt.xlim(mn,mx)
		plt.ylim(0, 0.04)
		plt.legend((p1,p2), ('Images with no timestamp in EXIF',
							 'Image timestamped ' + str(DD.date())),
				   prop=dict(size=9))
		for f in formats:
			plt.savefig(prefix + f)

	prefix = 'exif8'
	if 8 in makeplots:
		plt.clf()
		plt.subplots_adjust(**spargs)
		p2 = plt.semilogy(C.times, ptime, '0.5', lw=3)
		p1 = plt.semilogy(C.times, C.ptime, 'k-')
		print 'p1', p1
		print 'p2', p2
		p1 = p1[0]
		p2 = p2[0]

		ltimes = [datetime.datetime(x/12, x%12+1, 1)
				  for x in range(2007*12 + 5, 2008*12 + 4, 1)]
		plt.xticks([datetomjd(d) for d in ltimes],
				   [str(d.date()) if ((i%2) == 0) else ''
					for i,d in enumerate(ltimes)],
				   rotation=15, verticalalignment='top',
				   horizontalalignment='right')
		plt.xlabel('date')
		plt.ylabel('probability density (d$^{-1}$)')
		plt.xlim(mn,mx)
		plt.legend((p1,p2), ('no EXIF date; $p_{\mathrm{emp}}(t)$',
				     'EXIF ' + str(DD)),
			   'upper right', prop=dict(size=12))
		plt.ylim(6.e-5,1.5e0)
		for f in formats:
			plt.savefig(prefix + f)

	prefix = 'exif9'
	if 9 in makeplots:
		plt.clf()
		plt.subplots_adjust(**spargs)
		plt.hist(np.fmod([d for (wcs,d) in data if d > 0], 1.0),
				 range=(0, 1), bins=24, histtype='step', color='k')
		plt.xlabel('time of day (fractional MJD)')
		plt.ylabel('number of images')
		for f in formats:
			plt.savefig(prefix + f)


	xyz = C.xyz_at_times()
	I = []
	T = []
	wcs = Tan()

	nomjd = 0
	notin = 0

	for ii, (w, mjd) in enumerate(data):
		if mjd == 0:
			nomjd += 1
			continue
		wcs.set(*w)
		J = C.find_points_in_wcs(wcs, xyz)
		if len(J) == 0:
			notin += 1
			continue
		if not all(np.diff(J) == 1):
			print 'not contiguous'
			print J
		tlo,thi = min(C.times[J]), max(C.times[J])
		I.append(ii)
		T.append((tlo, thi, mjd))
	I = np.array(I)
	T = np.array(T)

	print 'No MJD:', nomjd
	print 'Not in:', notin
	print 'T:', len(I)

	prefix = 'exif10'
	if 10 in makeplots:
		plt.clf()
		plt.subplots_adjust(**spargs)

		Tlo = T[:,0]
		Thi = T[:,1]
		Texif = T[:,2]

		DT = (Thi - Tlo) / 2.
		print 'max comet-in-image span:', 2*DT.max()
		exifdt = Texif - (Thi + Tlo)/2.
		print 'EXIF diff range:', exifdt.min(), exifdt.max()

		def nlmap(X):
			return np.arcsinh(X)
		#np.sign(exifdt) * np.log10(np.abs(exifdt))

		J = np.argsort(DT)
		DT = DT[J]
		exifdt = exifdt[J]

		Xf = np.arange(len(J))
		# reverse index
		Ir = np.arange(len(J)-1, -1, -1)
		Xr = Ir

		# envelope
		Yf = nlmap(DT)
		env  = np.vstack((np.append(np.append([Xf[0]-10], Xf), [Xf[-1]+10]),
						  np.append(np.append([Yf[0]],    Yf), [Yf[-1]]))).T
		Yf = -nlmap(DT[Ir])
		renv  = np.vstack((np.append(np.append([Xr[0]+10], Xr), [Xr[-1]-10]),
						   np.append(np.append([Yf[0]],    Yf), [Yf[-1]]))).T
		plt.gca().add_artist(Polygon(np.vstack((env, renv)),
									 closed=True, fc='0.9', ec='0.7'))
		mn,mx = -8,8
		margin = (mx-mn)*0.005
		arrowlen = (mx-mn)*0.01
		K = np.flatnonzero(nlmap(exifdt) > mx)
		print len(K), 'above scale'
		if len(K):
			print 'exifdt', exifdt[K]
			#plt.plot(K, (mx)*np.ones_like(K), 'r.')
			for k in K:
				#plt.arrow(k, mx-(arrowlen+margin), 0, arrowlen, color='r')

				#arrowlen = (mx-mn)*0.05
				#plt.arrow(k, mx-margin, 0, -(arrowlen + margin), color='r',
				#		  head_starts_at_zero=False,
				#		  length_includes_head=True,
				#		  overhang=0, shape='right',
				#		  head_width=1, head_length=1)

				arrowlen = (mx-mn)*0.03
				plt.annotate('', (k,mx-margin), xytext=(k,mx-margin-arrowlen),
							 arrowprops=dict(color='k', lw=1,  headwidth=4, frac=0.5,
											 width=0.5))
								 
		K = np.flatnonzero(nlmap(exifdt) < mn)
		print len(K), 'below scale'
		if len(K):
			plt.plot(K, mn*np.ones_like(K), 'r.')

		y0 = nlmap(exifdt-0.5)
		y1 = nlmap(exifdt+0.5)
		# enforce minimum size in plot space
		minfrac = 0.007
		yminsize = (mx - mn)*minfrac
		I = np.abs(y1 - y0) < yminsize

		ym = (y0 + y1)/2.
		y0[I] = ym[I] - yminsize/2.
		y1[I] = ym[I] + yminsize/2.
		plt.plot(np.vstack((Xf,Xf)),
				 np.vstack((y0,y1)),
				 '-', color='0.2', lw=1.5)
				 #'k-', alpha=0.8, lw=1.5)

		#J = np.logical_not(I)
		#plt.plot(np.vstack((Xf[J],Xf[J])),
		#		 np.vstack((y0[J],y1[J])),
		#		 'k-', alpha=0.8, lw=1.5)
		#xminsize = len(J) * minfrac
		#plt.plot(np.vstack((Xf[I]-xminsize/2, Xf[I]+xminsize/2)),
		#		 np.vstack((ym[I],ym[I])),
		#		 'k-', alpha=0.8, lw=1.5)

		#plt.plot(nlmap(exifdt), 'k.', ms=1, alpha=1)

		plt.xlabel('image number (sorted by comet traversal duration)')
		plt.ylabel('EXIF time - Midpoint of comet-in-image time range (days)')
		yt = np.array([-1000, -100, -10, -1, 0, 1, 10, 100, 1000, 10000])
		plt.yticks(nlmap(yt.astype(float)), yt)
		plt.xlim(-0.5, len(J)-0.5)
		plt.ylim(mn, mx)
		for f in formats:
			plt.savefig(prefix + f)

		print 'exif10:', len(J)



def old_exifplots(data):

	if True:
		# EXIF dates plot.
		T = []
		xyz = C.xyz_at_times()
		di = 0
		print 'main(): n data', len(data)
		nodates = 0
		yesdates = 0
		notin = 0
		yesin = 0
		for ii, (w, mjd) in enumerate(data):
			if mjd == 0:
				nodates += 1
				continue
			yesdates += 1
			J = C.find_points_in_wcs(w, xyz)
			tlo,thi = 0,0
			if len(J) == 0:
				notin += 1
				c = 'r'
			else:
				yesin += 1
				tlo,thi = min(C.times[J]), max(C.times[J])
				c = 'b'
			T.append((w, mjd, tlo, thi))
			di += 1
		print 'main(): no dates', nodates
		print 'main(): yes dates', yesdates
		print 'main(): not in', notin
		print 'main(): yes in', yesin

		allsizes = np.array([w.pixel_scale()*np.sqrt(w.imagew*w.imageh) for w,d in data])
		I = np.argsort(allsizes)
		alldates = np.array([d for w,d in data])

		allsizes = allsizes[I]
		alldates = alldates[I]
		imgnum = np.arange(len(I))

		plt.clf()
		#plt.semilogx(allsizes[I]/3600., alldates[I], 'r.')
		plt.plot(imgnum, alldates, 'k.')
		J = np.flatnonzero(alldates == 0)
		#plt.plot(J, np.zeros_like(J) + 54200, 'r.')
		plt.ylabel('MJD')
		plt.xlabel('Image number, ordered by size')
		plt.xlim(0, len(allsizes))
		chunk = 50
		c0 = range(0, len(allsizes), chunk)
		fracexif = []
		for j in c0:
			dates = alldates[j:j+chunk]
			fracexif.append(float(np.sum(dates > 0)) / len(dates))
		fracexif = np.array(fracexif)
		#plt.twinx()
		cstep = np.array(zip(c0, c0[1:] + [len(allsizes)])).ravel()
		plt.plot(cstep, 54200. + fracexif.repeat(2) * 100, 'r-')
		plt.axhline(54300., color='r', linestyle='--')
		#plt.plot(I, allsizes[I]/3600., 'r.')
		#plt.ylabel('image size (deg)')
		plt.ylim(54200, 54700)
		plt.savefig('exif-0.png')

		#sys.exit(0)

		MJD = np.array([t[1] for t in T])
		Tlo = np.array([t[2] for t in T])
		Thi = np.array([t[3] for t in T])
		imsize = np.array([t[0].pixel_scale()*np.sqrt(t[0].imagew*t[0].imageh) for t in T])

		plt.clf()
		I = np.argsort(MJD)
		for i,(mjd,tlo,thi) in enumerate(zip(MJD[I], Tlo[I], Thi[I])):
			if tlo == 0:
				c = 'r'
			else:
				plt.plot([tlo, thi], [i, i], 'k-', alpha=0.25)
				c = 'k'
			plt.plot([mjd], [i], '.', color=c)
		plt.xlim(54200, 54700)
		plt.savefig('exif-1.png')

		plt.clf()
		I = np.argsort(imsize)
		for i,(mjd,tlo,thi) in enumerate(zip(MJD[I], Tlo[I], Thi[I])):
			if tlo == 0:
				c = 'r'
			else:
				plt.plot([tlo, thi], [i, i], 'k-', alpha=0.25)
				c = 'k'
			plt.plot([mjd]*2, [i-0.5,i+0.5], '-', color=c)
		plt.xlim(54200, 54700)
		plt.savefig('exif-2.png')

		I = np.argsort(imsize)
		I1 = I[:len(I)/3]
		I2 = I[len(I)/3:len(I)*2/3]
		I3 = I[len(I)*2/3:]
		plt.clf()
		for k,I in enumerate([I1,I2,I3]):
			plt.subplot(1, 3, k+1)
			J = np.argsort(MJD[I])
			J = I[J]
			for i,(mjd,tlo,thi) in enumerate(zip(MJD[J], Tlo[J], Thi[J])):
				if tlo == 0:
					c = 'r'
				else:
					plt.plot([tlo, thi], [i, i], 'k-', alpha=0.25)
					c = 'k'
				plt.plot([mjd]*2, [i-0.5,i+0.5], '-', color=c)
			plt.xlim(54200, 54700)
		plt.savefig('exif-3.png')

		plt.clf()
		I = np.argsort(imsize)
		J = (Tlo[I] > 0)
		I = I[J]
		for i,(mjd,tlo,thi) in enumerate(zip(MJD[I], Tlo[I], Thi[I])):
			tmid = (tlo + thi)/2.
			plt.plot([tlo-tmid, thi-tmid], [i, i], 'k-', alpha=0.25)
			plt.plot([mjd-tmid]*2, [i-0.5,i+0.5], '-', color=c)
		plt.xlim(-30, 30)
		plt.savefig('exif-4.png')

		plt.clf()
		I = np.argsort(imsize)
		J = (Tlo[I] > 0)
		I = I[J]
		for i,(mjd,tlo,thi,ii) in enumerate(zip(MJD[I], Tlo[I], Thi[I], I)):
			t0 = mjd
			plt.plot([tlo-t0, thi-t0], [i, i], 'k-', alpha=0.25)
			#plt.plot([mjd-t0, mjd-t0], [i-0.5, i+0.5], 'k-')
		plt.axvline(0, color='k')
		plt.xlim(-30, 30)
		plt.savefig('exif-5.png')

		plt.clf()
		plt.errorbar(MJD, (Tlo+Thi)/2., yerr=(Thi-Tlo)/2., fmt='.', ecolor='k')
		plt.xlabel('EXIF')
		plt.ylabel('Comet in image')
		plt.xlim(54200, 54700)
		plt.ylim(54200, 54700)
		plt.savefig('exif-6.png')


if __name__ == '__main__':
	if not os.path.exists(pfn):
		print 'Run mcmc_comet.py --exif-plots to regenerate pickle file', pfn
		sys.exit(0)

	print sys.argv
	args = sys.argv[1:]
	make = None
	if len(args):
		make = [int(a) for a in args]
		print 'Making plots', make

	D = unpickle_from_file(pfn)
	exifplots(*D, makeplots=make)
