import matplotlib
if __name__ == '__main__':
	matplotlib.use('Agg')
import os
import fnmatch
import sys
import pyfits

import numpy as np
import matplotlib
import pylab as plt
from math import floor
from glob import glob

from astrometry.util.util import *
from astrometry.util.pyfits_utils import *
from astrometry.blind.plotstuff import *
from astrometry.util.starutil_numpy import *

#   wget "http://www.astro.princeton.edu/~dstn/temp/scihackday/wcsfiles.tgz"
#   mkdir wcs-files
#   cd wcs-files/
#   tar xzf ../wcsfiles.tgz 
#   cd ..
#   wget "http://www.astro.princeton.edu/~dstn/temp/scihackday/wcsdata.txt"
#   text2fits.py wcsdata.txt wcsdata.fits

def make_pixdensity_map():
	T = fits_table('wcsdata.fits')
	T.pixdensity = 1./(T.pixscale**2)
	I = np.argsort(T.pixdensity)

	# pixels per square degree
	ppsd = (3600./T.pixscale)**2
	print 'pixdensity range', ppsd[I[0]], ppsd[I[-1]]
	#pixdensity range 11.0038206766 9650346614.7

	W = 2000
	H = W/2
	# Create Hammer-Aitoff WCS of the appropriate size.
	wcs = anwcs_create_allsky_hammer_aitoff(0., 0., W, H)
	plot = Plotstuff(outformat='png', size=(W, H))
	plot.wcs = wcs
	outline = plot.outline

	# Copied from holmes/wcsplots.py : pixdensity()
	fn1pat = 'pixdensity-%02i.png'
	ii = 0
	b = -1
	scales = []
	# This determines how much truncation error occurs
	# vs how many rounds are required.
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
			fn = os.path.join('wcs-files', T.wcsfn[i])
			print fn, 'ppsd', ppsd[i], '-> inc', inc, '-> sum', thissum
			if thissum > 1.:
				(mr,mg,mb,ma) = plotstuff_get_maximum_rgba(plot.pargs)
				print 'actual max red value:', mr
				thissum = float(mr) / 256. + inc
				print 'reset sum to', thissum
				if thissum > 1.:
					break
			ii += 1
			plot.rgb = (inc, 0, 0)
			plot_outline_set_wcs_file(outline, fn, 0)
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

	pyfits.writeto('psum.fits', psum)
	T = tabledata()
	T.scales = scales
	T.writeto('scales.fits')

	plt.clf()
	plt.imshow(psum, interpolation='nearest', origin='lower')
	plt.hot()
	plt.savefig('psum.png')




if __name__ == '__main__':
	if not os.path.exists('psum.fits'):
		make_pixdensity_map()

	psum = pyfits.open('psum.fits')[0].data
	H,W = psum.shape
	vmin,vmax = 2.5, 5.5
	I = (np.log10(psum+1.) - vmin) / (vmax - vmin)
	cmap = matplotlib.cm.hot
	rgb = cmap(I)
	rgb = np.clip(rgb * 255, 0, 255).astype(np.ubyte)

	wcs = anwcs_create_allsky_hammer_aitoff(0., 0., W, H)
	plot = Plotstuff(outformat='png', size=(W, H))
	plot.wcs = wcs
	img = plot.image
	img.set_image_from_numpy(rgb)
	plot.plot('image')

	# grid
	plot.fontsize = 12
	ras = [-180, -120, -60, 0, 60, 120, 180]
	decs = [-60, -30, 0, 30, 60]
	# dark gray
	plot.rgb = (0.3,0.3,0.3)
	plot.apply_settings()
	for ra in ras:
		plot.line_constant_ra(ra, -90, 90)
		plot.stroke()
	for dec in decs:
		plot.line_constant_dec(dec, -180, 180)
		plot.stroke()
	plot.color = 'gray'
	plot.apply_settings()
	for ra in ras:
		plot.move_to_radec(ra, 0)
		plot.text_radec(ra, 0, '%i'%((ra+360)%360))
		plot.stroke()
	for dec in decs:
		if dec != 0:
			plot.move_to_radec(0, dec)
			plot.text_radec(0, dec, '%+i'%dec)
			plot.stroke()

	plot.color = 'green'
	l = 0
	for b in range(360):
		r,d = ecliptictoradec(l,b) #lbtoradec(l,b)
		if b == 0:
			plot.move_to_radec(r,d)
		else:
			plot.line_to_radec(r,d)
	plot.stroke()

	plot.write('pixdensity.png')

	#plt.imshow(np.log10(1e3+I[0].data), vmax=5.5)
