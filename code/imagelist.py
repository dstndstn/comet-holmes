from astrometry.net import django_commandline

from astrometry.net.tile.models import *

import pickle

if __name__ == '__main__':

	outs = []
	imgs = MapImageSet.objects.get(id=2).images.all()
	for img in imgs:
		rd = img.job.get_radec_center()
		if rd is None:
			continue
		(ra,dec) = rd
		r = img.job.get_field_radius()
		t0 = img.exif_orig_date
		t1 = img.exif_date

		outs.append((ra,dec,r,t0,t1))

	print pickle.dumps(outs)
