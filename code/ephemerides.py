if __name__ == '__main__':
	f = open('ephemerides.txt')

	print '''
	function get_holmes_path() {
	poly = [];
	'''

	for line in f:
		if line.startswith('#'):
			continue
		words = line.split()
		ra = float(words[4])
		ra *= 15.
		dec = float(words[5])
		#print ra, dec

		print 'poly.push(new GLatLng(%g,%g));' % (dec, 360.-ra)
		
	print '''
	gpoly = new GPolyline(poly, "#00FF88", 2, 0.8);
	return gpoly;
	}'''
