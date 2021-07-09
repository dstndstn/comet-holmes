# Copyright 2011 Dustin Lang (Princeton) and David W. Hogg (NYU).
# All rights reserved.

import os
from glob import glob
import hashlib
from astrometry.util.file import *
import numpy as np
#import shutil

def main():
	import optparse
	import sys
	
	parser = optparse.OptionParser(usage='%prog <input-dir> <output-dir>')
	parser.add_option('-w', '--wcs', dest='wcs', default=False,
					  action='store_true',
					  help='Only consider .jpg files that have corresponding .wcs?')

	opt,args = parser.parse_args()
	if len(args) != 2:
		parser.print_help()
		sys.exit(-1)

	indir = args[0]
	outdir = args[1]

	if not os.path.exists(outdir):
		print 'Creating dir', outdir
		os.makedirs(outdir)

	if opt.wcs:
		wcsfns = glob(os.path.join(indir, 'holmes-*.wcs'))
		jpegfns = [fn.replace('.wcs','.jpg') for fn in wcsfns]
		print 'found', len(wcsfns), 'WCS filenames'
	else:
		jpegfns = glob(os.path.join(indir, 'holmes-*.jpg'))
		print 'found', len(jpegfns), 'JPEG filenames'

	pnmhashes = []

	hashes = []
	for i,jpegfn in enumerate(jpegfns):
		if not os.path.exists(jpegfn):
			raise RuntimeError('No such file: %s' % jpegfn)
		m = hashlib.md5()
		m.update(read_file(jpegfn))
		md5 = m.hexdigest()
		print md5, jpegfn
		hashes.append(md5)

		if False:
			pnmfn = '/tmp/pnm'
			cmd = 'jpegtopnm %s > %s' % (jpegfn, pnmfn)
			rtn = os.system(cmd)
			if not (os.WIFEXITED(rtn) and os.WEXITSTATUS(rtn) == 0):
				raise RuntimeError('jpegtopnm failed')
			m = hashlib.md5()
			m.update(read_file(pnmfn))
			md5 = m.hexdigest()
			pnmhashes.append(md5)
		

	hashes = np.array(hashes)
	#print hashes
	U,I = np.unique(hashes, return_index=True)
	print len(U), 'unique JPEGs'

	if False:
		pnmhashes = np.array(pnmhashes)
		U2,I2 = np.unique(pnmhashes, return_index=True)
		print len(U2), 'unique PNMs'

	for i in I:
		jpegfn = jpegfns[i]
		J = np.flatnonzero(hashes == hashes[i])
		JJ = (J != i)
		J = J[JJ]
		if len(J):
			#print 'file', jpegfn
			#print 'duplicates:', [jpegfns[j] for j in J]
			#print jpegfn, ' '.join([jpegfns[j] for j in J])
			dups = [jpegfn] + [jpegfns[j] for j in J]
			dups.sort()
			print ' '.join(dups)

		if opt.wcs:
			wcsfn = wcsfns[i]
			cmd = 'cp %s %s %s' % (wcsfn, jpegfn, outdir)
		else:
			cmd = 'cp %s %s' % (jpegfn, outdir)
		rtn = os.system(cmd)
		if not (os.WIFEXITED(rtn) and os.WEXITSTATUS(rtn) == 0):
			raise RuntimeError('copy failed')
		

if __name__ == '__main__':
	main()
