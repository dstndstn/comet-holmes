import sys
from astrometry.util.file import *

if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) != 1:
		print 'Usage: %s <short-name>' % (sys.argv[0])
		sys.exit(-1)

	name = args[0]

	(imgs,siteimgs) = unpickle_from_file('ysiteimages-%s.pickle' % name)

	# each of 'imgs' and 'siteimgs' contains tuples: (img url, w, h, filesize)

	print >>sys.stderr, '%i images' % len(imgs)
	print >>sys.stderr, '%i site images' % len(siteimgs)
	uimgs = list(set(imgs + siteimgs))
	print >>sys.stderr, '%i unique images' % len(uimgs)
	uimgs.sort()

	imgs = []
	for i,(url,w,h,filesize) in enumerate(uimgs):
		fn = '%s-%04i' % (name, i)
		imgs.append((url, w, h, filesize))
		print '-O %s %s' % (fn, url)

	pickle_to_file(imgs, 'yimageurlfns-%s.pickle' % name)
	
