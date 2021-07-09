import matplotlib
if __name__ == '__main__':
	matplotlib.use('Agg')
import os
import fnmatch
import sys

import numpy as np
import matplotlib
import pylab as plt
from math import floor

#from astrometry.util.sip import *
from astrometry.util.util import *
from astrometry.blind.plotstuff import *
from astrometry.util import jpl
from astrometry.util.file import *
from glob import glob

def plotit(width=0, gridstep=0, outfn=None,
		   dirfn = '2010-04-16-holmes-dedup',
		   ra=None, dec=None, rot=None, size=(1000,800),
		   wcsfn=None,
		   imagealpha = 0.02, resample=True,
		   wcsoutfn=None,
		   lastN=0, firstN=0,
		   format=PLOTSTUFF_FORMAT_PNG,
		   subpct = 0.25,
		   ):
	plot = Plotstuff()
	ispdf = (format == PLOTSTUFF_FORMAT_PDF)
	if ispdf:
		if saveimgsfn:
			plot.outfn = saveimgsfn
		elif outfn:
			plot.outfn = outfn
		elif ephoutfn:
			plot.outfn = ephoutfn
		assert(plot.outfn is not None)
		
	plot.size = size
	#print 'plotit(): wcsfn=',wcsfn
	if wcsfn is not None:
		#print 'setting plot WCS:', wcsfn
		plot.wcs_file = wcsfn
	else:
		plot.set_wcs_box(ra, dec, width)
	if rot is not None:
		plot.rotate_wcs(rot)
	plot.outformat = format
	plot.color = 'black'
	plot.plot('fill')

	fns = glob(os.path.join(dirfn, '*.wcs'))
	if lastN != 0:
		fns = fns[-lastN:]
	if firstN != 0:
		fns = fns[:firstN]

	pargs = plot.pargs
	outline = plot.outline
	img = plot.image
	img.alpha = imagealpha
	plot.alpha = 1.

	plot.pargs.op = CAIRO_OPERATOR_ADD
	for wcsfn in fns:
		jpegfn = wcsfn.replace('.wcs', '.jpg')
		if not os.path.exists(jpegfn):
			raise 'No such file ' + jpegfn
		print jpegfn
		if resample:
			img.resample = 1
		else:
			img.resample = 0
		plot_image_set_wcs(img, wcsfn, 0)
		plot_image_set_filename(img, jpegfn)

		if subpct > 0:
			rgb = plot_image_get_percentile(plot.pargs, img, subpct)
			print '  percentile', (100.*subpct), '=', rgb
			negrgb = [-x for x in rgb]
			plot_image_add_to_pixels(img, negrgb)
		plot.plot('image')

	plot.write(outfn)
	del plot


def weighted_coadd(size, comet='holmes'):
	if comet == 'holmes':
		eph = 'holmes-ephem.txt'
		plotwcsfn = 'holmes-coadd-1.wcs'
		dirfn = '2010-04-16-holmes-dedup'
		fnpat = 'holmes-wc-%i.png'
		gridstep = 5
	elif comet == 'hya':
		eph = 'hya-eph.txt'
		plotwcsfn = 'hya-1.wcs'
		dirfn = '2010-10-02-hyakutake'
		fnpat = 'hya-wc-%i.png'
		gridstep = 20
	else:
		raise RuntimeError('no such comet ' + comet)
		
	pfn = '%s-wco3.pickle' % comet
	if os.path.exists(pfn):
		(w,wimg,uimg,umimg) = unpickle_from_file(pfn)
		#W = wimg / w[:,:,np.newaxis]
		#print W.shape
		plt.figure(figsize=(10,10), dpi=100)
		plt.clf()
		plt.subplots_adjust(left=0, bottom=0, top=1, right=1)
		# imshow(np.clip(np.arcsinh(W / 10)/3., 0, 1))
		# log10 weights peak around 4.5-5, min ~ 4, max ~ 8.
		# clf(); hist(log10(w.ravel()), 100)
		W2 = (wimg + 0.*1e5) / (w[:,:,np.newaxis] + 1e5)
		# clf(); imshow(np.clip(np.arcsinh(W2 / 10)/3., 0, 1))
		# clf(); imshow(np.clip((W2+10) / 100., 0, 1))

		if comet == 'holmes':
			plt.imshow(np.clip((W2) / 75., 0, 1))
		elif comet == 'hya':
			plt.imshow(np.clip((W2+5) / 40., 0, 1))
		coaddimg = '%s-wco-1.png' % comet
		plt.savefig(coaddimg)

		for fn in ['%s-wco-2.png'%comet, '%s-wco-2.pdf'%comet]:
			plot = Plotstuff()
			plot.size = size
			plot.wcs_file = plotwcsfn
			if fn.endswith('.png'):
				plot.outformat = PLOTSTUFF_FORMAT_PNG
			elif fn.endswith('.pdf'):
				plot.outformat = PLOTSTUFF_FORMAT_PDF
			plot.outfn = fn
			img = plot.image
			img.set_file(coaddimg)
			plot.plot('image')

			set_grid_style(plot)
			plot.plot_grid(rastep=gridstep, decstep=gridstep)
			plot.alpha = 0.5
			grid = plot.grid

			if comet == 'holmes':
				grid.ralabeldir = DIRECTION_NEG
				grid.ralo = 40
				grid.rahi = 65
				plot.plot_grid(rastep=0, decstep=0,
							   ralabelstep=gridstep, declabelstep=gridstep)
			elif comet == 'hya':
				plot.plot_grid(rastep=0, decstep=0,
							   ralabelstep=gridstep, declabelstep=0)
				griddec = np.arange(0., 90.1, 20.)
				rafordeclabels=240
				for dec in griddec:
					plot_grid_add_label(plot.pargs, rafordeclabels, dec, dec)
			plot_ephem(plot, eph)
			plot.write(fn)
			print 'Wrote', fn

		return

	fns = glob(os.path.join(dirfn, '*.wcs'))
	plot = Plotstuff()
	plot.size = size
	plot.wcs_file = plotwcsfn
	plot.outformat = PLOTSTUFF_FORMAT_PNG
	outline = plot.outline
	plotscale = anwcs_pixel_scale(plot.wcs)
	pixscale = np.array([Tan(wcsfn).pixel_scale() for wcsfn in fns])
	# pixels per square degree
	ppsd = (3600./pixscale)**2
	I = np.argsort(ppsd)
	ii = 0

	weight = None
	wimg = None
	uimg = None
	umimg = None

	plot.op = CAIRO_OPERATOR_OVER
	img = plot.image
	img.resample = 1
	J = np.argsort(-ppsd)
	ploti = 0
	for ri,i in enumerate(J):
		plot.op = CAIRO_OPERATOR_OVER
		plot.color = 'black'
		plot.plot('fill')
		plot.op = CAIRO_OPERATOR_ADD

		scale = ppsd[i]
		wcsfn = fns[i]
		jpegfn = wcsfn.replace('.wcs', '.jpg')
		print ri, '/', len(J)
		sfac = plotscale / pixscale[i]
		print 'scale factor:', sfac
		if sfac > 2:
			sfac = int(floor(sfac))
			tmpfn = '/tmp/tmp.jpg'
			cmd = 'jpegtopnm %s | pnmscale -reduce %i | pnmtojpeg > %s' % (jpegfn, sfac, tmpfn)
			print cmd
			rtn = os.system(cmd)
			if not (os.WIFEXITED(rtn) and os.WEXITSTATUS(rtn) == 0):
				raise RuntimeError('jpegtopnm failed: rtn %i' % rtn)
			img.set_wcs_file(wcsfn, 0)
			anwcs_scale_wcs(img.wcs, 1./float(sfac))
			img.set_file(tmpfn)
		else:
			img.set_wcs_file(wcsfn, 0)
			img.set_file(jpegfn)
		#ok,ra,dec,r = anwcs_get_radec_center_and_radius(plot.image.wcs)
		#print 'center', ra, dec
		#ok,x,y = anwcs_radec2pixelxy(plot.wcs, ra, dec)
		#print 'x,y', x,y
		plot_image_read(plot.pargs, img)
		plot_image_add_to_pixels(img, [1,1,1])
		plot.plot('image')
		#plot.write('hwc-%04i.png' % ri)
		im = plot.get_image_as_numpy()
		#print 'im max', im.max()
		im = im[:,:,:3].astype(float)
		print 'im max', im.max()
		#print 'im shape', im.shape
		if im.max() == 0:
			print 'No pixels copied'
			continue

		# we added one to the pixel values to detect areas that were painted
		Wacc = scale * (im[:,:,0] > 0)
		# subtract it back out.
		im -= 1.

		imacc = np.empty_like(im)
		for c in range(3):
			imacc[:,:,c] = (im[:,:,c] - np.median(im[:,:,c])) * scale

		if wimg is None:
			wimg = imacc
		else:
			wimg += imacc

		if uimg is None:
			uimg = im
		else:
			uimg += im
		if umimg is None:
			umimg = imacc / scale
		else:
			umimg += imacc / scale

		if weight is None:
			weight = Wacc
		else:
			weight += Wacc

		#print 'Wacc range', W.min(), W.max()
		#print 'weight range', weight.min(), weight.max()

		if ri % 10 == 0:
			plt.clf()
			plt.imshow(wimg / wimg.max())
			plt.savefig('%s-wc-%03i-g.png' % (comet, ploti))

			plt.clf()
			plt.imshow(weight / weight.max())
			plt.savefig('%s-wc-%03i-h.png' % (comet, ploti))

			plt.clf()
			print 'max weight:', weight.max()
			WW = wimg / np.maximum(1, weight[:,:,np.newaxis])
			plt.imshow(WW / WW.max())
			plt.savefig('%s-wc-%03i-i.png' % (comet, ploti))
			ploti += 1

	pickle_to_file((weight, wimg, uimg, umimg), pfn)


def holmes_coadd_figs(ra, dec, width, size):
	eph = 'holmes-ephem.txt'
	plotfn = 'holmes-coadd-1.png'
	plotwcsfn = 'holmes-coadd-1.wcs'
	
	if not os.path.exists(plotfn):
		# alpha = 0.03, subpct = 0.3
		plotit(ra=ra, dec=dec, width=width, size=size, outfn=plotfn,
			   gridstep=1, resample=True,
			   imagealpha = 0.03, subpct = 0.3)

	if not os.path.exists(plotwcsfn):
		plot = Plotstuff()
		plot.size = size
		plot.set_wcs_box(ra, dec, width)
		anwcs_write(plot.wcs, plotwcsfn)

	fullsz = (1000,1000)
	halfsz = (500,500)

	coaddimg = plotfn

	halfimg = 'holmes-coadd-1b.png'
	cmd = 'pngtopnm %s | pnmscale -reduce 2 | pnmtopng > %s' % (coaddimg, halfimg)
	print 'Running:', cmd
	os.system(cmd)

	for fn,fullsize in [('holmes-coadd-2.png',True), ('holmes-coadd-2.pdf',True),
						('holmes-coadd-2b.pdf',False),
						]:
		plot = Plotstuff()
		if fullsize:
			plot.size = fullsz
		else:
			plot.size = halfsz
		plot.wcs_file = plotwcsfn
		if not fullsize:
			anwcs_scale_wcs(plot.wcs, 0.5)
			
		if fn.endswith('.png'):
			plot.outformat = PLOTSTUFF_FORMAT_PNG
		elif fn.endswith('.pdf'):
			plot.outformat = PLOTSTUFF_FORMAT_PDF
		plot.outfn = fn
		img = plot.image
		if fullsize:
			img.set_file(coaddimg)
		else:
			img.set_file(halfimg)
		plot.plot('image')

		set_grid_style(plot, not fullsize)
		gridstep = 5
		plot.plot_grid(rastep=gridstep, decstep=gridstep)
		plot.alpha = 0.5
		grid = plot.grid
		grid.ralabeldir = DIRECTION_NEG
		grid.ralo = 40
		grid.rahi = 65
		plot.plot_grid(rastep=0, decstep=0,
					   ralabelstep=gridstep, declabelstep=gridstep)
		plot_ephem(plot, eph, not fullsize)

		plot.write(fn)
		print 'Wrote', fn

def sumpngs(fns):
	Isum = None
	for fn in fns:
		I = plt.imread(fn)
		# cut from rgba to r.
		I = I[:,:,0]
		#print 'unique vals:', np.around(255. * np.unique(I))
		if Isum is None:
			Isum = I
		else:
			Isum += I
	print 'Sum:', Isum.shape, Isum.min(), Isum.max()
	# convert back to integers
	Isum = (Isum * 255.).astype(int)
	return Isum


def holmes_ndensity_figs(ra, dec, w, size):
	# number of sub-plots
	B = 6
	# If I were being more cunning I would pack the filled WCSes and
	# the edges into two planes of the PNG.
	fn1pat = 'holmes-footprints-1-%i.png'
	fn2pat = 'holmes-footprints-2-%i.png'
	fn1s = [fn1pat % b for b in range(B)]
	fn2s = [fn2pat % b for b in range(B)]
	if not all([os.path.exists(fn) for fn in fn1s+fn2s]):
		plot = Plotstuff()
		plot.size = (1000,1000)
		plot.set_wcs_box(ra, dec, w)
		plot.outformat = PLOTSTUFF_FORMAT_PNG
		eta = 1./256.
		outline = plot.outline
		fns = glob('2010-04-16-holmes-dedup/*.wcs')
		# do them in batches of Nb to prevent saturation.
		Nb = 250
		b = 0
		while len(fns):
			for fill,pat in [(1, fn1pat), (0, fn2pat)]:
				plot.op = CAIRO_OPERATOR_OVER
				plot.color = 'black'
				plot.plot('fill')
				plot.op = CAIRO_OPERATOR_ADD
				plot.rgb = (eta, eta, eta)
				for wcsfn in fns[:Nb]:
					print wcsfn
					plot_outline_set_wcs_file(outline, wcsfn, 0)
					outline.fill = fill
					plot.plot('outline')
				plot.write(pat % b)
			b += 1
			fns = fns[Nb:]

	setup_plot()
	Isum = sumpngs(fn1s)
	outsum = sumpngs(fn2s)
	print 'Max number of images:', Isum.max()
	plt.imshow(Isum)
	plt.hot()
	#plt.colorbar(cax=colorbar_axes(plt.gca(), frac=0.08))
	gridra,griddec = get_holmes_grid()
	add_grid(ra, dec, w, size, 'holmes-footprints-1.wcs',
			 gridra, griddec, gridlines=False)
	plt.savefig('holmes-footprints-1.pdf')
	clim = plt.gci().get_clim()

	plt.clf()
	plt.imshow(Isum + 5*outsum)
	plt.clim(*clim)
	plt.hot()
	eph = 'holmes-ephem.txt'
	add_grid(ra, dec, w, size, 'holmes-footprints-1.wcs',
			 gridra, griddec, ephem=eph)
	plt.savefig('holmes-footprints-2.pdf')


def numberdensity(fns, fn1pat, rdw, size, wcsfn=None):
	plot = Plotstuff()
	plot.size = (1000,1000)
	if rdw is not None:
		plot.set_wcs_box(ra, dec, w)
	elif wcsfn is not None:
		plot.wcs_file = wcsfn
	else:
		# need WCS!
		assert(False)
	# do them in batches of Nb to prevent saturation.
	Nb = 250
	plot.outformat = PLOTSTUFF_FORMAT_PNG
	outline = plot.outline
	outline.fill = 1
	eta = 1./256.
	b = 0
	while len(fns):
		print
		print 'starting round', b
		plot.op = CAIRO_OPERATOR_OVER
		plot.color = 'black'
		plot.plot('fill')
		plot.op = CAIRO_OPERATOR_ADD
		plot.rgb = (eta, 0, 0)
		for wcsfn in fns[:Nb]:
			print wcsfn
			plot_outline_set_wcs_file(outline, wcsfn, 0)
			plot.plot('outline')
		plot.write(fn1pat % b)
		b += 1
		fns = fns[Nb:]
	B = b
	# read the images and sum them
	psum = None
	for b in range(B):
		fn = fn1pat % b
		I = plt.imread(fn)
		# cut from rgba to r.
		I = I[:,:,0].astype(float) * 255.
		if psum is None:
			psum = I
		else:
			psum += I
	return psum

def pixdensity(fns, fn1pat, rdw, size, wcsfn=None):
	plot = Plotstuff()
	plot.size = (1000,1000)
	if rdw is not None:
		plot.set_wcs_box(ra, dec, w)
	elif wcsfn is not None:
		plot.wcs_file = wcsfn
	else:
		# need WCS!
		assert(False)
	plot.outformat = PLOTSTUFF_FORMAT_PNG
	outline = plot.outline

	pixscale = []
	for wcsfn in fns:
		wcs = Tan(wcsfn)
		pixscale.append(wcs.pixel_scale())
	pixscale = np.array(pixscale)
	# pixels per square degree
	ppsd = (3600./pixscale)**2
	I = np.argsort(ppsd)
	ii = 0
	b = -1
	scales = []
	# This determines how much truncation error occurs
	# vs how many rounds (fn1pat pngs) are required.
	eta = 5./256.

	while True:
		if ii >= len(I):
			break
		b += 1
		print
		print 'starting round', b
		thissum = 0.
		i = I[ii]
		scale = ppsd[i]
		print 'min pixels per square degree:', scale
		scales.append(scale)
		plot.op = CAIRO_OPERATOR_OVER
		plot.color = 'black'
		plot.plot('fill')
		plot.op = CAIRO_OPERATOR_ADD
		outline.fill = 1
		while ii < len(I):
			i = I[ii]
			inc = ppsd[i] / scale * eta
			thissum += inc
			print fns[i], 'ppsd', ppsd[i], '-> inc', inc, '-> sum', thissum
			if thissum > 1.:
				(mr,mg,mb,ma) = plotstuff_get_maximum_rgba(plot.pargs)
				print 'actual max red value:', mr
				thissum = float(mr) / 256. + inc
				print 'reset sum to', thissum
				if thissum > 1.:
					break
			ii += 1
			plot.rgb = (inc, 0, 0)
			plot_outline_set_wcs_file(outline, fns[i], 0)
			plot.plot('outline')
		plot.write(fn1pat % b)

	scales = np.array(scales)

	psum = None
	for b,S in enumerate(scales):
		fn = fn1pat % b
		I = plt.imread(fn)
		# cut from rgba to r.
		I = I[:,:,0].astype(float) * S
		if psum is None:
			psum = I
		else:
			psum += I
	return psum


def set_grid_style(plot, halfsize=False):
	plot.text_bg_alpha = 0
	plot.rgb = (0., 1., 0.)
	plot.alpha = 0.3
	#plot.lw = 0.5
	if halfsize:
		plot.lw = 0.5
		plot.fontsize = 10
	else:
		plot.lw = 1.
		plot.fontsize = 20

def set_ephem_style(plot, halfsize=False):
	plot.rgb = (0., 0.7, 0.)
	plot.alpha = 1.0
	if halfsize:
		plot.lw = 1.5
	else:
		plot.lw = 3.

def plot_ephem(plot, eph, halfsize=False):
	(ra,dec,jd) = jpl.parse_radec(open(eph).read())
	set_ephem_style(plot, halfsize)
	plotstuff_move_to_radec(plot.pargs, ra[0], dec[0])
	for r,d in zip(ra,dec):
		plotstuff_line_to_radec(plot.pargs, r, d)
	plotstuff_stroke(plot.pargs)

def holmes_pixdensity_figs(ra, dec, w, size):
	# Pixels per square degree: max 4178620.97762
	# min 255.227250345
	picfn = 'pixdensity.pickle'
	fns = glob('2010-04-16-holmes-dedup/*.wcs')
	if not os.path.exists(picfn):
		fn1pat = 'holmes-footprints-3-%i.png'
		size = (1000,1000)
		psum = pixdensity(fns, fn1pat, (ra,dec,w), size)
		pickle_to_file(psum, picfn)
	else:
		psum = unpickle_from_file(picfn)

	print 'Pixels per square degree: max', psum.max()
	print 'min', psum.min()

	vmin,vmax = 2.5, 6.5
	I = (np.log10(psum) - vmin) / (vmax - vmin)
	cmap = matplotlib.cm.hot
	rgb = cmap(I)
	rgb = np.clip(rgb * 255, 0, 255).astype(np.ubyte)

	matplotlib.rc('font', family='computer modern roman')
	matplotlib.rc('font', size=14)
	matplotlib.rc('text', usetex=True)
	plt.figure(figsize=(0.5,4))
	plt.clf()
	ax = plt.gca()
	print 'ax pos', ax.get_position()
	#ax.set_position([0, 0, 1, 1])
	ax.set_position([0.01, 0., 0.35, 1.])
	#cax,kw = matplotlib.colorbar.make_axes(ax, fraction=0.9)
	#fraction 0.15; fraction of original axes to use for colorbar
	#pad 0.05 if vertical, 0.15 if horizontal; fraction of original axes between colorbar and new image axes
	#shrink 1.0; fraction by which to shrink the colorbar
	#aspect 20; ratio of long to short dimensions

	cb = matplotlib.colorbar.ColorbarBase(
		ax, cmap=cmap,
		norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax),
		)
	cb.set_ticks([3,4,5,6], update_ticks=False)
	cb.set_ticklabels(['$10^3$', '$10^4$', '$10^5$', '$10^6$'], update_ticks=False)
	cb.update_ticks()
	plt.savefig('cbar.png')
	plt.savefig('cbar.pdf')

	#values=None, boundaries=None, orientation='vertical', extend='neither', spacing='uniform', ticks=None, format=None, drawedges=False, filled=True)

	plotwcsfn = 'holmes-footprints-1.wcs'
	eph = 'holmes-ephem.txt'

	for fn,grid,ann in [('holmes-footprint-1.png', False, False),
						('holmes-footprint-1.pdf', False, False),
						('holmes-footprint-2.png', True,  True),
						('holmes-footprint-2.pdf', True,  True),
						('holmes-footprint-3.png', True,  False),
						('holmes-footprint-3.pdf', True,  False),]:
		plot = Plotstuff()
		plot.size = (1000,1000)
		plot.wcs_file = plotwcsfn
		if fn.endswith('.png'):
			plot.outformat = PLOTSTUFF_FORMAT_PNG
		elif fn.endswith('.pdf'):
			plot.outformat = PLOTSTUFF_FORMAT_PDF
		plot.outfn = fn
		img = plot.image
		img.set_image_from_numpy(rgb)
		plot.plot('image')

		if grid:
			set_grid_style(plot)
			gridstep = 5
			plot.plot_grid(rastep=gridstep, decstep=gridstep)
			plot.alpha = 0.5
			grid = plot.grid
			grid.ralabeldir = DIRECTION_NEG
			grid.ralo = 40
			grid.rahi = 65
			plot.plot_grid(rastep=0, decstep=0,
						   ralabelstep=gridstep, declabelstep=gridstep)
		if ann:
			plot_ephem(plot, eph)

		plot.write(fn)
		print 'Wrote', fn

def hya_figs():
	eph = 'hya-eph.txt'
	plotwcsfn = 'hya-1.wcs'
	halfwcsfn = 'hya-1-half.wcs'
	dirnm = '2010-10-02-hyakutake'

	ra,dec = 180.,55.
	width = 140.
	W,H = 1000,1000
	size = (W,H)
	#pixscale = width / W

	gridra,griddec = np.arange(0., 360.1, 20.), np.arange(0., 90.1, 20.)
	rafordeclabels=240

	import pyfits
	for w,h,fn in [(W,H,plotwcsfn), (W/2, H/2, halfwcsfn)]:
		pixscale = width / w
		hdu = pyfits.PrimaryHDU(None)
		hdr = hdu.header
		hdr.update('CTYPE1', 'RA---ARC')
		hdr.update('CTYPE2', 'DEC--ARC')
		hdr.update('WCSAXES', 2)
		hdr.update('EQUINOX', 2000.)
		hdr.update('LONPOLE', 180.)
		hdr.update('LATPOLE', 0.)
		hdr.update('CRVAL1', ra)
		hdr.update('CRVAL2', dec)
		hdr.update('CRPIX1', w/2.)
		hdr.update('CRPIX2', h/2.)
		hdr.update('CUNIT1', 'deg')
		hdr.update('CUNIT2', 'deg')
		hdr.update('CD1_1', pixscale)
		hdr.update('CD1_2', 0.)
		hdr.update('CD2_2', pixscale)
		hdr.update('CD2_1', 0.)
		hdr.update('IMAGEW', w)
		hdr.update('IMAGEH', h)
		hdu.writeto(fn, clobber=True)

	# Just the images
	plotfn = 'hya-coadd-1.png'
	# alpha=0.02, subpct=0 -- not quite saturated; quite washed out.
	# alpha=0.03, subpct=0.1 -- hya-coadd-3a.png -- nice; bit dark
	# alpha=0.07, subpct=0.1 -- just right
	if not os.path.exists(plotfn):
		plotit(outfn=plotfn, size=size, resample=True,
			   wcsfn=plotwcsfn, dirfn=dirnm, format=PLOTSTUFF_FORMAT_PNG,
			   imagealpha=0.07, subpct = 0.1)

	coaddimg = plotfn

	# halfsize image
	halfimg = 'hya-coadd-1b.png'
	cmd = 'pngtopnm %s | pnmscale -reduce 2 | pnmtopng > %s' % (coaddimg, halfimg)
	print 'Running:', cmd
	os.system(cmd)

	fullsz = (1000,1000)
	halfsz = (500,500)

	for fn,fullsize in [('hya-coadd-2.png',True), ('hya-coadd-2.pdf',True),
						('hya-coadd-2b.pdf', False)]:
		plot = Plotstuff()
		if fullsize:
			plot.size = fullsz
			plot.wcs_file = plotwcsfn
		else:
			plot.size = halfsz
			plot.wcs_file = halfwcsfn
		if fn.endswith('.png'):
			plot.outformat = PLOTSTUFF_FORMAT_PNG
		elif fn.endswith('.pdf'):
			plot.outformat = PLOTSTUFF_FORMAT_PDF
		plot.outfn = fn
		img = plot.image
		if fullsize:
			img.set_file(coaddimg)
		else:
			img.set_file(halfimg)
		plot.plot('image')

		set_grid_style(plot, not fullsize)
		gridstep = 20
		plot.plot_grid(rastep=gridstep, decstep=gridstep)
		plot.alpha = 0.5
		grid = plot.grid
		plot.plot_grid(rastep=0, decstep=0,
					   ralabelstep=gridstep, declabelstep=0)
		for dec in griddec:
			plot_grid_add_label(plot.pargs, rafordeclabels, dec, dec, '%.0f')
		plot_ephem(plot, eph, not fullsize)

		plot.write(fn)
		print 'Wrote', fn
		

	# Pixel density map
	size = (1000,1000)
	picfn = 'hya-pixdensity.pickle'
	fns = glob(os.path.join(dirnm, '*.wcs'))
	if not os.path.exists(picfn):
		fn1pat = 'hya-footprints-3-%i.png'
		psum = pixdensity(fns, fn1pat, None, size, wcsfn=plotwcsfn)
		pickle_to_file(psum, picfn)
	else:
		psum = unpickle_from_file(picfn)

	if False:
		# Maximum number density: 152.
		fnpat = 'hya-footprints-4-%i.png'
		nsum = numberdensity(fns, fnpat, None, size, wcsfn='hya-coadd-1.wcs')
		print 'Maximum number density:', nsum.max()

	vmin,vmax = 0., 4.
	I = (np.log10(np.maximum(1, psum)) - vmin) / (vmax - vmin)
	cmap = matplotlib.cm.hot
	rgb = cmap(I)
	rgb = np.clip(rgb * 255, 0, 255).astype(np.ubyte)

	for fn,ann in [('hya-footprint-1.png', False),
				   ('hya-footprint-1.pdf', False),
				   ('hya-footprint-2.png', True),
				   ('hya-footprint-2.pdf', True)]:
		plot = Plotstuff()
		plot.size = (1000,1000)
		plot.wcs_file = plotwcsfn
		if fn.endswith('.png'):
			plot.outformat = PLOTSTUFF_FORMAT_PNG
		elif fn.endswith('.pdf'):
			plot.outformat = PLOTSTUFF_FORMAT_PDF
		plot.outfn = fn
		img = plot.image
		img.set_image_from_numpy(rgb)
		plot.plot('image')

		if ann:
			set_grid_style(plot)
			gridstep = 20
			plot.plot_grid(rastep=gridstep, decstep=gridstep)
			plot.alpha = 0.5
			grid = plot.grid
			plot.plot_grid(rastep=0, decstep=0,
						   ralabelstep=gridstep, declabelstep=0)
			for dec in griddec:
				plot_grid_add_label(plot.pargs, rafordeclabels, dec, dec, '%.0f')
			plot_ephem(plot, eph)

		plot.write(fn)
		print 'Wrote', fn



if __name__ == '__main__':
	log_init(2)
	fits_use_error_system()

	# Holmes fig
	ra,dec = 52.5, 45
	# 35 deg wide in RA
	w = 35 * np.cos(np.deg2rad(dec))
	size = (1000, 1000)

	holmes_coadd_figs(ra, dec, w, size)


	sys.exit(0)
	hya_figs()


	holmes_pixdensity_figs(ra, dec, w, size)

	sys.exit(0)


	weighted_coadd(size, comet='holmes')
	weighted_coadd(size, comet='hya')
	


	#holmes_ndensity_figs(ra, dec, w, size)



	# how many images are outside this plot area?
	dirfn = '2010-04-16-holmes-dedup'
	fns = glob(os.path.join(dirfn, '*.wcs'))
	wcs1 = anwcs_create_box(ra, dec, w, size[0], size[1])
	nin = 0
	nout = 0
	for fn in fns:
		print fn
		wcs2 = anwcs_open(fn, 0)
		if anwcs_overlaps(wcs1, wcs2, 100):
			nin += 1
		else:
			nout += 1
	print 'N overlapping:', nin
	print 'N outside:', nout
	sys.exit(0)
