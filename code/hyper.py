import matplotlib
matplotlib.use('Agg')
import numpy as np
import pylab as plt
import os
from astrometry.util.file import *
from astrometry.util.plotutils import antigray

from astrometry.util import celestial_mechanics as cm
from mcmc_comet import *


def xax(ns=2000, firsttick=500, step=500):
	plt.xlabel('MCMC iteration')
	plt.xticks(np.arange(firsttick,ns+1, step))
	plt.xlim(-10,ns+10)

def plotparam(x, A, eps=False, sub=5):
	if eps:
		plt.plot(x, A, '-', color='0.8', lw=0.5, zorder=-10)
	else:
		plt.plot(x, A, 'k-', alpha=0.05) #, zorder=10)
	y = np.mean(A, axis=1)
	yerr=np.std(A, axis=1)
	x = np.append(x[::sub], x[-1])
	y = np.append(y[::sub], y[-1])
	yerr = np.append(yerr[::sub], yerr[-1])
	#print 'x', x
	#print 'y', y
	#print 'yerr', yerr
	plt.errorbar(x, y, yerr, color='k', barsabove=True) #, zorder=20)

def gauss(X, mu, sigma):
	#return 1./(2.*pi*sigma**2) * np.exp(e)
	# prevent underflow in exp
	e = -(X - mu)**2 / (2.*sigma**2)
	rtn = np.zeros_like(e)
	I = (e > -700)  # -> ~ 1e-304
	rtn[I] = 1./(2.*pi*sigma**2) * np.exp(e[I])
	return rtn

# 1-D
def em_step(X, mu, sigma, A):
	# E:
	# fg = p( Y, Z=f | theta ) = p( Y | Z=f, theta ) p( Z=f | theta )

	N = len(X)
	K = len(A)
	D = 1
	print 'N', N, 'D', D, 'K', K

	fg = np.zeros([N, K])
	for i in range(len(A)):
		print 'mu', mu[i], 'sig', sigma[i]
		fg[:,i] = gauss(X, mu[i], sigma[i]) * A[i]

	# normalize:
	# fore = p( Z=f | Y, theta )
	fore = fg / np.sum(fg, axis=1)[:,np.newaxis]
	print 'fg', fg.shape
	print 'fore', fore.shape
	print fore.sum(axis=1)
	assert(all(np.isfinite(fore.ravel())))

	# M:
	# maximize mu, sigma:
	for i in range(len(A)):
		mu[i] = np.sum(fore[:,i] * X) / np.sum(fore[:,i])

		# 2.*sum(fore) because X,mu are 2-dimensional.
		sigma[i] = np.sqrt(np.sum(fore[:,i] * (X - mu[i])**2) / np.sum(fore[:,i]))

	# maximize A.
	A = np.mean(fore, axis=0)
	print 'new A', A

	# avoid multiplying 0 * -inf = NaN
	#I = (fg > 0)
	#lfg = np.zeros_like(fg)
	#lfg[I] = np.log(fg[I])
	# Total expected log-likelihood
	#Q = np.sum(fore*lfg + back*lbg)
	return (mu, sigma, A)


def hyperparam_plots(lnp, params):
	print 'Hyperparam plots'
	pgood = params[:,:,6]
	pexif = params[:,:,7]
	lncent = params[:,:,8]

	(ns,nw) = pgood.shape
	print 'steps', ns, 'walkers', nw

	#mx = np.mean(lnp, axis=1)[-1]
	mx = np.max(lnp)
	dlnp = lnp - mx

	HYPS = [
		(pgood, 'pgood', r'$p_{\mathrm{good}}$', None),
		(pexif, 'pexif', r'$p_\mathrm{EXIF}$', None),
		(np.exp(-lncent), 'eta', r'$\eta$', (0.3,0.5)),
		#(np.exp(-lncent), 'eta', r'$\eta$', (0.33,0.40)),
		(dlnp, 'lnp', 'ln posterior probability (relative)', (-50,10)),
		#(dlnp, 'lnp', 'ln posterior probability (relative)', (-20,1)),
		]

	xaxp = dict(ns=ns, firsttick=1000, step=1000)
	
	for i,(param,base,label,ylim) in enumerate(HYPS):

		print 'Hyperparam', base

		sub=(ns / 20)
		plt.clf()
		plt.subplots_adjust(**spargs)
		plotparam(step, param, sub=sub)
		xax(**xaxp)
		plt.ylabel(label)
		if ylim is not None:
			plt.ylim(*ylim)
		plt.savefig('hyp-' + base + '.png')
		plt.savefig('hyp-' + base + '.pdf')

		plt.clf()
		plt.subplots_adjust(**spargs)
		plotparam(step, param, eps=True, sub=sub)
		xax(**xaxp)
		if ylim is not None:
			plt.ylim(*ylim)
		plt.ylabel(label)
		plt.savefig('hyp-' + base + '.eps')

	for i,(param,base,label,ylim) in enumerate(HYPS):
		for j,(jparam,jbase,jlabel,jylim) in enumerate(HYPS):
			if j <= i:
				continue
			print 'Hyperparams', base, jbase

			sub=100
			plt.clf()
			#plt.plot(param.ravel(), jparam.ravel(), '.', alpha=0.02)
			hargs = {}
			if ylim is not None and jylim is not None:
				hargs.update(range=(ylim, jylim))
			(H,xe,ye) = np.histogram2d(param.ravel(), jparam.ravel(), 100,
									   **hargs)
			ima = dict(extent=(min(xe), max(xe), min(ye), max(ye)),
					   aspect='auto', cmap=antigray,
					   interpolation='nearest', origin='lower')
			plt.imshow(H.T, **ima)
			plt.gray()
			plt.subplots_adjust(**spargs)
			plt.xlabel(label)
			plt.ylabel(jlabel)
			if ylim is not None:
				plt.xlim(*ylim)
			if jylim is not None:
				plt.ylim(*jylim)
			plt.savefig('hh-%s-%s.png' % (base,jbase))



def xyz_plots(params, epoch, color=True, usealpha=True):
	print 'xyz_plots'
	# Plot in x,y,z space
	times = np.linspace(epoch-400, epoch+400, 200)

	ecliptic=False

	C = CometMCMCxv()
	C.set_true_params()
	txyzs = C.heliocentric_xyz_at_times(times=times, ecliptic=ecliptic)
	print 'xyzs', txyzs.shape

	mxyz = C.heliocentric_xyz_at_times(times=[epoch], ecliptic=ecliptic)
	xc,yc,zc = mxyz[0,:]

	ldays = [(2007, 8, 1), (2007, 10, 1), (2007, 12, 1),
			 (2008, 2, 1), (2008, 4, 1),]
	ltimes = [datetomjd(datetime.datetime(*x)) for x in ldays]
	lxyz = C.heliocentric_xyz_at_times(times=ltimes, ecliptic=ecliptic)

	exyz = EMB_xyz_at_times(times)
	elxyz = EMB_xyz_at_times(ltimes)

	ns,nw,nparm = params.shape

	print 'Number of samples:', nw
	# CUT
	nw = 32

	sxyzs = []
	for i in range(nw):
		C.set_params(params[-1,i,:])
		C.set_epoch(epoch)
		sxyzs.append(C.heliocentric_xyz_at_times(times=times,
												 ecliptic=ecliptic))
	# R=1.5
	R=2.

	#spargs2 = dict(bottom=0.15, left=0.15, right=0.95, top=0.95)
	spargs2 = dict(bottom=0.11, left=0.11, right=0.99, top=0.99)
	for fig in range(3):
		plt.figure(num=fig+1, figsize=(4,4), dpi=dpi)
		plt.subplots_adjust(**spargs2)
		plt.clf()

	labels = [r'$\boldsymbol{%s}$ (AU)'%s for s in 'xyz']
	lims = [(xc-R, xc+R), (yc-R, yc+R), (zc-R, zc+R)]

	for i,xyzs in enumerate([txyzs, exyz] + sxyzs):

		# Comet path: samples
		if usealpha:
			stys = [dict(color='k', alpha=0.1, zorder=10)]
		else:
			stys = [dict(color='0.4', zorder=10, lw=0.5)]

		if i == 0:
			# True comet path
			if color:
				stys = [dict(color='r', lw=2, alpha=1., zorder=12)]
				lsty = dict(color='r', mec='r', mfc='r', ms=5, zorder=12)
			else:
				if usealpha:
					stys = [dict(color='k', lw=1, alpha=0.5, zorder=11),
							dict(color='k', lw=2, ls='--', alpha=1., zorder=12)]
				else:
					stys = [dict(color='0.5', lw=1, zorder=11),
							dict(color='k', lw=2, ls='--', alpha=1., zorder=12)]
				lsty = dict(color='k', mec='k', mfc='k', ms=5, zorder=12)
		elif i == 1:
			# EMB path
			if usealpha:
				stys = [dict(color='k', alpha=0.25, lw=2, zorder=10)]
			else:
				stys = [dict(color='0.75', lw=2, zorder=10)]
			if usealpha:
				elsty = dict(color='k', mec='k', mfc='k', ms=5, zorder=12, alpha=0.5)
			else:
				elsty = dict(color='0.5', mec='k', mfc='k', ms=5, zorder=12)


			

		for j,(xi,yi) in enumerate(([0,2], [0,1], [2,1])):
			plt.figure(num=j+1)
			for sty in stys:
				plt.plot(xyzs[:,xi], xyzs[:,yi], '-', **sty)

			dx,dy = R/20., R/40.
			if j == 2:
				dy *= -1
			if j == 1:
				dx *= 1.5

			if color:
				tc = 'r'
			else:
				tc = 'k'

			if i == 0:
				plt.plot(lxyz[:,xi], lxyz[:,yi], 'o', **lsty)
				for x,y,d in zip(lxyz[:,xi], lxyz[:,yi], ldays):
					alignargs = dict(ha='left', va='center')
					ang = 0.
					if j == 0:
						alignargs = dict(ha='left', va='bottom', rotation=45.)
					plt.text(x+dx, y+dy, '%i-%02i-%02i' % d, 
							 fontdict=dict(size=10, color=tc), **alignargs)
					#bbox=dict(facecolor='w', alpha=1.0, edgecolor='w', fill=True, zorder=11))
					#bbox=dict(boxstyle='square', facecolor='0.9', alpha=1.0, edgecolor='k', fill=True, zorder=11))
				plt.xlabel(labels[xi])
				plt.ylabel(labels[yi])
				plt.xlim(lims[xi])
				plt.ylim(lims[yi])

			elif i == 1:
				plt.plot(elxyz[:,xi], elxyz[:,yi], 'o', **elsty)
				if usealpha:
					plt.plot(np.vstack((elxyz[:,xi], lxyz[:,xi])),
							 np.vstack((elxyz[:,yi], lxyz[:,yi])), '-', color='k', lw=1, alpha=0.25)
				else:
					plt.plot(np.vstack((elxyz[:,xi], lxyz[:,xi])),
							 np.vstack((elxyz[:,yi], lxyz[:,yi])), '-', color='0.75', lw=1)

	for fig in range(3):
		plt.figure(num=fig+1)
		
		ax = plt.axis()
		plt.xticks([0,1,2,3])
		plt.yticks([0,1,2,3])
		plt.axis(ax)
		for suff in ['png', 'pdf', 'eps']:
			if color:
				bwstr = ''
			else:
				bwstr = '-bw'
			plt.savefig('xyzorbit-%i%s.%s' % (fig, bwstr, suff))


def param_plots(params, epoch):
	print 'param_plots'

	print 'Params', params.shape
	# cut
	#nkeep = 800
	#sub = 1
	#params = params[-nkeep::sub, :, :]
	#print '->    ', params.shape

	spargs = dict(bottom=0.12, left=0.1, right=0.95, top=0.95)
	plt.subplots_adjust(**spargs)

	hyps = params[:,:,6:]
	# -lneta -> eta
	hyps[:,:,2] = 1./np.exp(hyps[:,:,2])

	for i,pname in enumerate([
		r'$\boldsymbol{p_{\mathrm{good}}}$',
		r'$\boldsymbol{p_{\mathrm{EXIF}}}$',
		r'$\boldsymbol{\eta}$',]):

		param = hyps[:,:,i]
		parm = param.ravel()
		mn = np.mean(parm)
		sd = np.std(parm)
		plt.clf()
		n,b,p = plt.hist(parm, 50, histtype='step', color='k', zorder=20)
		hint = sum(n)*(b[1]-b[0])
		plt.xlabel(pname)

		if i == 2:
			plt.xlim(0.34, 0.40)

		ax = plt.axis()
		ym = ax[3]
		d2 = '$%.2f \pm %.2f$'
		d3 = '$%.3f \pm %.3f$'
		fmts = [d3, d2, d3]

		yt,yl = plt.yticks()
		plt.yticks(yt, ['']*len(yt))
		plt.ylabel('fraction of samples')

		y1 = 0.2 * ym
		ax = plt.axis()
		x = np.linspace(ax[0], ax[1], 300)

		if i == 2:

			print 'parm', parm.shape
			X = parm
			#mu = np.array([2.6, 2.75, 2.85])
			#sigma = np.array([0.05, 0.05, 0.05])
			mu = np.array([0.35, 0.365, 0.385])
			sigma = np.array([0.005, 0.005, 0.005])
			A = np.array([0.33, 0.34, 0.33])

			for step in range(20):
				#print 'mu, sigma, A', mu, sigma, A
				print 'mu', mu
				print 'sigma', sigma
				print 'A', A
				mu,sigma,A = em_step(X, mu, sigma, A)

			y = np.zeros_like(x)
			for mui,si,ai in zip(mu,sigma,A):
				y += ai * 1./(np.sqrt(2.*np.pi)*si)*np.exp(-0.5*(mui-x)**2/si**2)
			#plt.plot(x, y*hint, 'k-', alpha=0.5, lw=2, zorder=15)
			plt.plot(x, y*hint, '-', color='0.5', lw=2, zorder=15)

			for j,(mui,si,ai) in enumerate(zip(mu,sigma,A)):
				y1 = 0.15 * ym
				if j == 1:
					y1 = 0.25 * ym
				plt.text(mui, y1, fmts[i] % (mui, si),
						 va='center', ha='center', color='k', zorder=25,
						 bbox=dict(facecolor='w', alpha=1.0, edgecolor='w', fill=True, zorder=24))


		else:
			plt.text(mn, y1, fmts[i] % (mn, sd),
					 va='center', ha='center', color='k')
			y = 1./(np.sqrt(2.*np.pi)*sd)*np.exp(-0.5*(mn-x)**2/sd**2)
			plt.plot(x, y*hint, '-', color='0.5', lw=2, zorder=15) #alpha=0.5, lw=2, zorder=15)

		if i == 0:
			ax = plt.axis()
			plt.xticks([0.86, 0.88, 0.90, 0.92, 0.94])
			plt.axis(ax)

		plt.savefig('hparam-hist-%i.png' % i)
		plt.savefig('hparam-hist-%i.pdf' % i)
		plt.savefig('hparam-hist-%i.eps' % i)

	print 'done hparams'

	ns,nw,nparm = params.shape
	X = params[:,:,:3]
	V = params[:,:,3:6]
	elements = np.zeros((ns, nw, 6))

	print 'converting to orbital elements...'
	for s in range(ns):
		for w in range(nw):
			elements[s,w,:] = cm.orbital_elements_from_phase_space_coordinates(X[s,w,:], V[s,w,:], cm.GM_sun)
	print 'done'

	a = elements[:,:,0]
	e = elements[:,:,1]
	I = elements[:,:,2]
	Om = elements[:,:,3]
	pom = elements[:,:,4]
	M = elements[:,:,5]

	C = CometMCMCxv()
	C.set_true_params()
	etrue = C.get_elements_at_time(epoch)
	#print 'True elements:', etrue

	PNAMES = [r'semi-major axis, $\boldsymbol{a}$ (AU)',
			  r'eccentricity, $\boldsymbol{e}$',
			  r'inclination, $\boldsymbol{I}$ (rad)',
			  r'long.~of asc.~node, $\boldsymbol{\Omega}$ (rad)',
			  r'argument of periapsis, $\boldsymbol{\omega}$ (rad)',
			  r'mean anomaly, $\boldsymbol{M}$ (rad)']

	xaxp = dict(firsttick=1000, step=1000)

	for i,pname, in enumerate(PNAMES):

		param = elements[:,:,i]

		if False:
			plt.clf()
			sub = (ns/20)
			plotparam(step, param, sub=sub)
			xax(**xaxp)
			plt.ylabel(pname)
			plt.savefig('param%i.png' % i)

		parm = param.ravel()
		mn = np.mean(parm)
		sd = np.std(parm)

		print 'mean', mn, 'std', sd
		print 'true', etrue[i]

		plt.clf()
		n,b,p = plt.hist(parm, 50, histtype='step', color='k', zorder=20)
		hint = sum(n)*(b[1]-b[0])
		plt.axvline(etrue[i], color='r', ls='--', lw=2) #color='r', alpha=0.5, lw=2)
		#ax = plt.axis()
		#x = np.linspace(b[0], b[-1], 300)
		#x = np.linspace(ax[0], ax[1], 300)
		#y = 1./(np.sqrt(2.*np.pi)*sd)*np.exp(-0.5*(mn-x)**2/sd**2)
		#plt.plot(x, y*hint, 'k-', alpha=0.5, lw=2, zorder=15)
		plt.xlabel(pname)

		ax = plt.axis()
		rx = ax[1]-ax[0]
		ym = ax[3]
		d2 = '$%.2f \pm %.2f$'
		d3 = '$%.3f \pm %.3f$'
		fmts = [d2,d2, d3, '$%.2f \pm %.3f$' ,d2,d2 ]
		d2 = '$%.2f$'
		d3 = '$%.3f$'
		tfmts = [d2,d2,d3, d2,d2,d2 ]
		thas = ['left', 'right', 'left', 'right', 'left', 'right']

		y1 = 0.2 * ym
		if i == 2:
			y1 = 0.2 * ym
		plt.text(mn, y1, fmts[i] % (mn, sd),
				 va='center', ha='center', color='k')
		dx = 0.02 * rx
		if thas[i] == 'right':
			dx *= -1
		if i == 2:
			y1 = 0.3 * ym
			dx = 0.01 * rx
		plt.text(etrue[i] + dx, y1, tfmts[i] % etrue[i],
				 va='center', ha=thas[i], color='r')

		y2 = 0.1 * ym
		dy2 = 0.02 * ym
		plt.plot([mn,mn],[y2-dy2,y2+dy2], '-', color='0.5', lw=2, alpha=1., zorder=10)
		plt.plot([mn,etrue[i]],[y2,y2], '-', color='0.5', lw=2, alpha=1., zorder=10)
		nsig = np.abs(mn-etrue[i]) / sd
		print 'nsigma', nsig

		lfracs = [0.5, 0.4, 0.4, 0.7, 0.7, 0.5]

		plt.text(mn + (etrue[i]-mn)*lfracs[i], y2, '$%.1f \sigma$' % nsig, color='k',
				 va='center', ha='center', zorder=13,
				 bbox=dict(facecolor='w', alpha=1.0, edgecolor='w', fill=True, zorder=12))

		#y3 = 0.05 * ym
		#plt.arrow(mn, y3, etrue[i]-mn, 0, label='$%.1f \sigma$' % nsig, color='k')

		yt,yl = plt.yticks()
		plt.yticks(yt, ['']*len(yt))
		plt.ylabel('fraction of samples')

		if i == 1:
			# eccentricity
			plt.xticks([0.30, 0.35, 0.40, 0.45])
			plt.xlim(0.27, 0.45)
		elif i == 2:
			# inclination
			plt.xticks([0.335, 0.340, 0.345])
			plt.xlim(0.332, 0.3453)
		elif i == 3:
			# long of asc node
			plt.xlim(5.59, 5.71)
		elif i == 4:
			# arg of peri
			plt.xticks([0.4, 0.5, 0.6, 0.7, 0.8])
			plt.xlim(0.39, 0.81)
		elif i == 5:
			# mean anomaly
			plt.xlim(0.45, 0.57)

		ax = plt.axis()
		#x = np.linspace(b[0], b[-1], 300)
		x = np.linspace(ax[0], ax[1], 300)
		y = 1./(np.sqrt(2.*np.pi)*sd) * np.exp(-0.5*(mn-x)**2/sd**2)
		plt.plot(x, y*hint, '-', color='0.5', lw=2, zorder=15) #'k-', alpha=0.5, lw=2, zorder=15)

		plt.savefig('param-hist-%i.png' % i)
		plt.savefig('param-hist-%i.pdf' % i)
		plt.savefig('param-hist-%i.eps' % i)

	from astrometry.util.plotutils import antigray

	
	ns,nw,nel = elements.shape
	nsam = ns*nw
	P = np.zeros((nel, nsam))
	for i in range(nel):
		P[i,:] = elements[:,:,i].ravel()
	print 'P:', P.shape
	C = np.cov(P)
	print 'C:', C
	Cinv = np.linalg.inv(C)
	mu = np.mean(P, axis=1)
	print 'mu', mu
	print mu.shape
	print 'etrue:', etrue
	print '  ', np.array(etrue)
	print '  ', np.array(etrue).shape
	d = mu - np.array(etrue)[:len(mu)]
	print 'd', d
	Cinvd = np.dot(Cinv, d)
	print 'Cinvd', Cinvd
	M = np.dot(d, Cinvd)
	print 'M', M

	V = np.var(P, axis=1)
	print 'V', V
	print '1/V', 1./V
	print 'chi2 = ', np.sum(d**2/V)

	spargs = dict(bottom=0.12, left=0.2, right=0.95, top=0.95)
	plt.subplots_adjust(**spargs)

	for i,pname, in enumerate(PNAMES):
		for j,pname2, in enumerate(PNAMES):
			if j <= i:
				continue
			parami = elements[:,:,i].ravel()
			paramj = elements[:,:,j].ravel()

			ti,tj = etrue[i], etrue[j]
			rng = ((min(ti, min(parami)), max(ti, max(parami))),
				   (min(tj, min(paramj)), max(tj, max(paramj))))
			H,xe,ye = np.histogram2d(parami, paramj, bins=100, range=rng)
			plt.clf()
			myargs = dict(extent=(min(xe), max(xe), min(ye), max(ye)),
						  aspect='auto',
						  interpolation='nearest', origin='lower')
			plt.imshow(H.T, cmap=antigray, **myargs)

			plt.plot(etrue[i], etrue[j], 'ro')
			plt.xlabel(pname)
			plt.ylabel(pname2)
			
			plt.savefig('params-%i-%i.png' % (i,j))
			print 'covariance', i,j




if __name__ == '__main__':
	fnpat = 'mcmc-default-nw064-bw008-emp-%04i.pickle'
	dirnm = '.'
	fnpat = os.path.join(dirnm, fnpat)

	epoch = 54416.0
	print 'Epoch', epoch
	print mjdtodate(epoch)
	
	params = []
	lnp = []
	step = []

	burn = 1000
	# steps = np.arange(burn, 2000+1, 10)
	# steps = np.arange(burn, 5000+1, 1)
	#steps = np.arange(burn, 5000+1, 1)
	#steps = np.arange(0, 5000+1, 1)
	steps = np.arange(0, 4000+1, 1)

	for k in steps:
		fn = fnpat % k
		print fn
		if not os.path.exists(fn):
			break
		(k,pos,lnprob,state) = unpickle_from_file(fn)
		params.append(pos)
		lnp.append(lnprob)
		step.append(k)
	
	step = np.array(step)
	params = np.array(params)
	lnp = np.array(lnp)
	print 'params', params.shape
	# (nsteps, nwalkers, nparams)
	print 'lnp', lnp.shape
	
	# dpi=300
	dpi=100
	plt.figure(figsize=(4,4), dpi=dpi)
	spargs = dict(bottom=0.1, left=0.2, right=0.95, top=0.9)
	matplotlib.rc('font', family='computer modern roman')
	matplotlib.rc('font', size=14)
	matplotlib.rc('text', usetex=True)
	matplotlib.rc('text.latex', preamble=r'\usepackage{bm}')
	matplotlib.rc('savefig', dpi=dpi)

	# hh-, hyp-
	#hyperparam_plots(lnp, params)

	# hparam-hist-*, param-hist-*, params-*-*.png
	#param_plots(params[burn:,:,:], epoch)

	# xyzorbit-*
	# For AJ .eps figs:
	xyz_plots(params[burn:,:,:], epoch, usealpha=False)
	#xyz_plots(params[burn:,:,:], epoch, color=False)

	sys.exit(0)


