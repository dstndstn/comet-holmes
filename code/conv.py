import os
for c in ['g','h','i']:
	j = 0
	os.system('rm %s???.jpg' % c)
	#for i in range(0, 102):
	for i in range(0, 74):
		fn = 'holmes-wc-%03i-%s.png' % (i,c)
		if not os.path.exists(fn):
			continue
		print fn
		os.system('pngtopnm %s | pnmtojpeg > %s%03i.jpg' % (fn, c, j))
		j += 1
	os.system('ffmpeg -i %s%%03d.jpg -r 10 -b 1000000 -y %s.mp4' % (c,c))
	os.system('rm %s???.jpg' % c)
