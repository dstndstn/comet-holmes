import sys

from pylab import *
from numpy import *
from scipy.ndimage.filters import *

img = imread('yimage-0980.png')
img = img.sum(axis=2)

clf()
imshow(img)
gray()
colorbar()
savefig('img.png')

clf()
hist(img.ravel(), 32)
savefig('hist.png')

#clf()
#blur = gaussian_filter(img, 10.)
#imshow(blur)
#savefig('blur.png')

clf()
med = median_filter(img, size=20.)
imshow(med)
savefig('median.png')

clf()
hist(med.ravel(), 32)
savefig('med-hist.png')

(h,w) = img.shape
gg = mgrid[0:h,0:w]
yy = gg[0]
xx = gg[1]

s = med.copy().ravel()
s.sort()

#clf()
#imshow(med)
#a=axis()
#for pct in arange(0.9, 1.0, 0.01):
#	thresh = s[min(len(s)-1, int(pct * len(s)))]
#	xt = median(xx[med >= thresh])
#	yt = median(yy[med >= thresh])
#	plot([xt],[yt], 'r.')
#	#text(xt, yt, '%.2f' % pct)
#thresh = s[min(len(s)-1, int(0.95 * len(s)))]

clf()
grad = gaussian_gradient_magnitude(med, 1.)
imshow(grad)
colorbar()
savefig('grad.png')

bestrxy = None
bestscore = -1e100

for r in range(20, 101, 10):
	s = r*2+4
	c = s/2
	gg = mgrid[0:s,0:s]
	yf = gg[0]
	xf = gg[1]
	filt = (abs(sqrt((yf-c)**2 + (xf-c)**2) - r) < 2.)
	clf()
	imshow(filt, interpolation='nearest')
	colorbar()
	savefig('filt-%03i.png' % r)
	hough = convolve(grad, filt)
	clf()
	imshow(hough)
	colorbar()
	savefig('hough-%03i.png' % r)
	print 'r=',r,' max hough=', hough.max()

	if hough.max() > bestscore:
		i = hough.argmax()
		xmax = xx.ravel()[i]
		ymax = yy.ravel()[i]
		bestrxy = (r, xmax, ymax)
		bestscore = hough.max()

(bestr,xc,yc) = bestrxy
print 'best r:', bestr

clf()
imshow(med)
a=axis()
ang=arange(0, 7, 0.1)
plot(sin(ang)*bestr + xc, cos(ang)*bestr + yc, 'r-')
axis(a)
savefig('best.png')
