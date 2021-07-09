from celestial_mechanics import *
from comet_likelihood import *
from astrometry.util import jpl
from astrometry.util.file import *
from numpy import *
from pylab import *

if __name__ == '__main__':

	jd,E = jpl.parse_orbital_elements('''2454101.500000000 = A.D. 2007-Jan-01 00:00:00.0000 (CT)
 EC= 1.670361927937051E-02 QR= 9.832911829245575E-01 IN= 9.028642170169823E-04
 OM= 1.762399911457168E+02 W = 2.867172565215373E+02 Tp=  2454104.323526433203
 N = 9.856169820212966E-01 MA= 3.572170843984495E+02 TA= 3.571221738148128E+02
 A = 9.999947139070439E-01 AD= 1.016698244889530E+00 PR= 3.652534468934519E+02''',
								   needSystemGM=False)
	GM = 2.9591310798672560E-04 #AU^3/d^2
	print 'Elements:', E
	(a,e,I,Omega,pomega,M0,nil) = E[0]

	dts = arange(365.) + 200
	dMdt = 	sqrt(GM / a**3)

	obsxyz = []
	for dt in dts:
		M = M0 + dMdt * dt
		(x,v) = phase_space_coordinates_from_orbital_elements(a,e,I,Omega,pomega,M,GM)
		obsxyz.append(x)
	obsxyz = array(obsxyz)
	print 'observer position', obsxyz.shape
	t0 = jd[0]
	times = t0 + dts

	#  comet_lnlikelihood(a, e, I, Omega, pomega, M0, t0
	#       times(K), observer(Kx3), imagexyz(Nx3), imagerad(N) )
	#

	ras,decs = meshgrid(arange(30,80,2), arange(30,60,2))

	#imgxyz = array([[1.,0,0]])
	#imgrad = array([1.])

	#print a,e,I,Omega,pomega,M, t0, times, obsxyz,imgxyz, imgrad

	jdc,Ec = jpl.parse_orbital_elements('''2454101.500000000 = A.D. 2007-Jan-01 00:00:00.0000 (CT)
 EC= 4.324200089863557E-01 QR= 2.053195945095214E+00 IN= 1.911310712121985E+01
 OM= 3.268678025165774E+02 W = 2.425019787429147E+01 Tp=  2454224.979476255365
 N = 1.432514716543015E-01 MA= 3.423113833073105E+02 TA= 3.137076507500682E+02
 A = 3.617456530538365E+00 AD= 5.181717115981516E+00 PR= 2.513063187712040E+03''',
								   needSystemGM=False)
	(a,e,I,Omega,pomega,M0,nil) = Ec[0]

	#x = comet_lnlikelihood(a,e,I,Omega,pomega,M, t0, times, obsxyz,
	#					   imgxyz, imgrad)
	lnls = []
	for ra,dec in zip(ras.ravel(), decs.ravel()):
		print 'ra,dec', ra,dec
		imgxyz = radectoxyz(ra, dec)
		imgrad = sqrt(deg2distsq(array([1.])))
		#print 'imgxyz', imgxyz.shape
		lnl = comet_lnlikelihood(a,e,I,Omega,pomega,M0, t0, times, obsxyz,
								 imgxyz, imgrad)
		lnls.append(lnl)

	lnls = array(lnls).reshape(ras.shape)
	print 'lnls:', lnls.min(), lnls.max()
	clf()
	for ra,dec,lnl in zip(ras.ravel(), decs.ravel(),lnls.ravel()):
		plot([ra], [dec], 'o', color=(exp(lnl - lnls.max()),0,0))

	(cras,cdecs,cjds) = jpl.parse_radec(read_file('../ephemeris/py/test-holmes-6.txt'))
	plot(cras, cdecs, 'r-')
	axis([30,80, 30,60])
	savefig('lnls.png')
	print x


