import sys
from urlparse import urlparse
from yahoo.search.image import ImageSearch
from util import *
from astrometry.util.file import *

if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) != 2:
		print 'Usage: %s <short-name> <query>' % (sys.argv[0])
		sys.exit(-1)

	name = args[0]
	query = args[1]


	# FIXME -- for site searches, we should probably be looking at
	# RefererUrl, if it exists, rather than the image url directly.


	appid = 'astrometry.net comet holmes project'
	#'lGICw2TV34HlhGrELJJdiQgo9E.SBko2PaUHIffvhZ6voJ8Vx9paBysWSi_YjyYHLFXt'
	#'Hx3duRXV34Fim3PxNexFBMQqVmWDnNzN23y4dneebBaKQBCRj10gip8_VFmfyI_dXXoHhKM_wOfq'
	start = 1
	nres = 50

	imgs = []
	while True:
		if start + nres > 1001:
			break
		# format="bmp","gif","jpeg","png"
		# site="www.yahoo.com"
		print 'Querying for:', query, 'starting with result', start
		search = ImageSearch(query=query, type='phrase',
							 format='jpeg',
							 app_id=appid, results=nres,
							 start=start, output='json')
		stream = search.open()
		json = stream.read()
		res = json2python(json)
		#niceprint(res)
		#sys.exit(0)
		
		rs = res['ResultSet']
		available = int(rs['totalResultsAvailable'])
		results = rs['Result']
		print '%i results available' % available
		for i,r in enumerate(results):
			url = r['Url']
			w = int(r['Width'])
			h = int(r['Height'])
			filesize = int(r['FileSize'])
			refurl = r['RefererUrl']
			imgs.append((url, w, h, filesize,refurl))
		start += len(results)
		if start >= available:
			break

	pickle_to_file(imgs, 'yimages-%s.pickle' % name)

	# Search within each site...

	sites = []
	for (url,w,h,filesize,refurl) in imgs:
		up = urlparse(url)
		sites.append(up.hostname)
	sites = list(set(sites))
	sites.sort()
	#print 'Sites:', '\n', '\n'.join(sites)

	siteimgs = []
	for si,s in enumerate(sites):
		start = 1
		nres = 50
		while True:
			if start + nres > 1001:
				break
			search = ImageSearch(query=query, type='phrase',
								 site=s, format='jpeg',
								 app_id=appid, results=nres,
								 start=start, output='json')
			print 'Searching site', s, si, 'of', len(sites)
			stream = search.open()
			json = stream.read()
			res = json2python(json)
			rs = res['ResultSet']
			available = int(rs['totalResultsAvailable'])
			results = rs['Result']
			print '%i results available' % available
			for i,r in enumerate(results):
				url = r['Url']
				w = int(r['Width'])
				h = int(r['Height'])
				filesize = int(r['FileSize'])
				refurl = r['RefererUrl']
				siteimgs.append((url, w, h, filesize,refurl))
			start += len(results)
			if start >= available:
				break

	pickle_to_file((imgs, siteimgs), 'ysiteimages-%s.pickle' % name)


