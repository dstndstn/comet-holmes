import matplotlib
matplotlib.use('Agg')
import pylab as plt

from mcmc_comet import *
from astrometry.util import celestial_mechanics as cm

def norm(x):
	return np.sqrt(np.sum(x**2, axis=1))

def main():
	C = CometMCMCxv()
	C.set_true_params()

	# days
	dt = 1.
	#C.set_times(datetomjd(datetime.datetime(2007, 3, 1)),
	#			datetomjd(datetime.datetime(2008, 9, 1)),
	#			dt)

	# tmin,tmax in the paper.
	C.set_times(datetomjd(datetime.datetime(2007, 7, 1)),
				datetomjd(datetime.datetime(2008, 5, 1)),
				dt)

	#C.set_times(datetomjd(datetime.datetime(2008, 4, 28)),
	#			datetomjd(datetime.datetime(2008, 5, 1)),
	#			dt)

	print 'c.x,v', C.x, C.v

	r1,d1 = C.radec_at_times(C.times, light_travel=True)
	r2,d2 = C.radec_at_times(C.times, light_travel=False)
	plt.clf()
	plt.plot(r1, d1, 'r-')
	plt.plot(r2, d2, 'b--')
	plt.savefig('ltt-radec.png')
	
	plt.clf()
	print r1.min(), r1.max(), d1.min(), d1.max()
	p1 = plt.plot(C.times, (r1-r2)*3600.*np.cos(np.deg2rad(d1)), 'r-')
	p2 = plt.plot(C.times, (d1-d2)*3600., 'b-')
	plt.legend((p1,p2), ('dRA*cos(Dec)', 'dDec'))
	plt.savefig('ltt-dradec.png')


	# the speed of light = 173.144483 Astronomical Units per day
	lightspeed = 173.144483

	exyz = C.earthxyz

	# test -- ~ one full period of the comet:
	t0 = datetomjd(datetime.datetime(2007, 7, 1))
	mm = []
	xyz1 = []
	TT = np.arange(t0, t0 + 6.88 * days_per_yr * 0.99, 10)
	for t in TT:
		E = C.get_elements_at_time(t)
		mm.append(E[5])
		x,dx = cm.orbital_elements_to_ss_xyz(E, light_travel=False)
		xyz1.append(x)
	print 'mm', max(mm), min(mm), max(mm)-min(mm)
	xyz1 = np.array(xyz1)
	plt.clf()
	plt.plot(xyz1[:,0], xyz1[:,1], 'r-')
	plt.savefig('comet.png')

	cxyz = []
	ltxyz = []
	ltxyz2 = []
	for e,t in zip(exyz, C.times):
		elements = C.get_elements_at_time(t)

		# (   a [AU]                   e                       i [rad]
		# [3.6184488483006279, 0.43256985249231866, 0.33357730745046937,
		#     Omega [rad]           pomega [rad]            M [rad]
		# 5.7048560402053807, 0.4236138975017858, 1.2107886403014576,
		#       GM [AU^3/day^2] )
		# 0.000295913107986725]
		
		x,dx = cm.orbital_elements_to_ss_xyz(elements, e, light_travel=False)
		cxyz.append(x)
		lx,ldx = cm.orbital_elements_to_ss_xyz(elements, e, light_travel=True)
		ltxyz.append(lx)

		# compute a rough light-travel correction by asking for the comet's position
		# at the initial light-travel-time ago.
		# [AU]
		dist0 = np.sqrt(np.sum(dx**2))
		print 'dist0:', dist0, 'AU'
		# [days]
		dt0 = dist0 / lightspeed
		print 'dt', dt0, 'days'
		print '  ', dt0*24.*60., 'minutes'
		el2 = C.get_elements_at_time(t - dt0)
		x2,dx2 = cm.orbital_elements_to_ss_xyz(el2, e, light_travel=False)
		ltxyz2.append(x2)

	cxyz = np.array(cxyz)
	ltxyz = np.array(ltxyz)
	ltxyz2 = np.array(ltxyz2)

	#print 'cxyz', cxyz.shape
	#print 'exyz', exyz.shape
	dx = cxyz - exyz
	ldx = ltxyz - exyz
	ldx2 = ltxyz2 - exyz
	#print dx.shape
	#print ldx.shape
	normdx = norm(dx)
	#print 'normdx', normdx.shape
	mind = np.min(normdx)
	print 'min dist:', mind, 'AU'
	print 'min light travel time, days:', mind / lightspeed

	#radius of Earth = 4.26349283x10-5 Astronomical Units
	re = 4.26e-5
	print 'max parallax due to Earth center vs surface:', np.rad2deg(re / mind), 'deg'

	dot = np.sum(dx * ldx, axis=1) / norm(dx) / norm(ldx)
	#print 'dot', dot.shape
	angle = np.rad2deg(np.arccos(dot))
	#print 'angle', angle.shape
	print 'max angle between no-LT and LT:', np.max(angle), 'deg'
	ltangle = angle
	dist = norm(dx - ldx)
	#print 'dist', dist.shape
	print 'max dist between no-LT and LT:', np.max(dist), 'AU'

	dot = np.sum(dx * ldx2, axis=1) / norm(dx) / norm(ldx2)
	angle = np.rad2deg(np.arccos(dot))
	print 'max angle between no-LT and LT-2:', np.max(angle), 'deg'
	ltangle3 = angle
	
	xyz0 = C.xyz_at_times(light_travel=False)
	xyz1 = C.xyz_at_times(light_travel=True)
	dot = np.sum(xyz0 * xyz1, axis=1) / norm(xyz0) / norm(xyz1)
	angle = np.rad2deg(np.arccos(dot))
	ltangle2 = angle
	print 'max angle between no-LT and LT:', np.max(angle), 'deg'

	ndx = dx / norm(dx)[:,np.newaxis]
	dx1 = ndx[1:,:]
	dx0 = ndx[:-1,:]
	dot = np.sum(dx1 * dx0, axis=1)
	angle = np.rad2deg(np.arccos(dot)) / dt
	print 'max angular speed:', np.max(angle), 'deg / day'
	angspeed = angle

	vc = np.diff(cxyz, axis=0)/dt
	speed = norm(vc)
	print 'max speed:', np.max(speed), 'AU/day'

	plt.clf()
	plt.plot(exyz[:,0], exyz[:,1], 'b-')
	plt.plot(cxyz[:,0], cxyz[:,1], 'r-')
	plt.plot(ltxyz[:,0], ltxyz[:,1], 'm-')
	plt.plot(ltxyz2[:,0], ltxyz2[:,1], 'g-')
	I = np.arange(0, cxyz.shape[0], 30)
	plt.plot(np.vstack((exyz[I,0], ltxyz[I,0])),
			 np.vstack((exyz[I,1], ltxyz[I,1])), 'm-')
	plt.plot(np.vstack((exyz[I,0], cxyz[I,0])),
			 np.vstack((exyz[I,1], cxyz[I,1])), 'k-')
	plt.axis('scaled')
	plt.savefig('par1.png')

	# angspeed: [deg/day]
	# aph: ["/hr]
	aph = angspeed * 3600. / 24.
	plt.clf()
	p1 = plt.plot(C.times[:-1], aph, 'r-')
	plt.ylabel('angular speed (arcsec/hr)')
	plt.xlabel('day (mjd)')
	plt.savefig('par3.png')

	ltthrs = normdx/lightspeed * 24.
	plt.clf()
	plt.plot(C.times, ltthrs, 'r-')
	plt.ylabel('light travel time (hrs)')
	plt.xlabel('day (mjd)')
	plt.savefig('par4.png')

	plt.clf()
	p1 = plt.plot(C.times, np.rad2deg(re / normdx) * 3600., 'r-')
	p2 = plt.plot(C.times, ltangle * 3600., 'b-')
	#p3 = plt.plot(C.times, ltangle2 * 3600., 'g-')
	#p4 = plt.plot(C.times, ltangle3 * 3600., 'm-')
	p5 = plt.plot(C.times[:-1], ltthrs[:-1] * aph, 'k-')
	plt.ylabel('size of effect (arcsec)')
	plt.xlabel('day (mjd)')
	plt.legend((p1,p2,p5), ('parallax', 'light travel', 'LTT * angular speed'))
	plt.savefig('par2.png')
	
	
if __name__ == '__main__':
	main()
