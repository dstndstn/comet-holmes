import math
import datetime
import pickle
import sys

from math import atan2, acos, sqrt, atan

from matplotlib.patches import Circle, Rectangle

import numpy.random

from pylab import *

import gettext
_= gettext.gettext

### Stolen from nightsky.py

piover180= math.pi/180.0
piunder180= 180.0/math.pi

def rad2deg(r):
	return r * piunder180

def deg2rad(d):
	return d * piover180

# class for handling all simple 3-d vector operations
class Vector:
	def __init__(self,start):
		self.data= start
	def __repr__(self):
		return self.__str__()
	def __str__(self):
		return '(%g, %g, %g)' % (self.data[0], self.data[1], self.data[2])
	def __neg__(self):
		return Vector([-self.data[0],-self.data[1],-self.data[2]])
	def __add__(self,other):
		return Vector([self.data[0]+other.data[0],
					   self.data[1]+other.data[1],
					   self.data[2]+other.data[2]])
	def __sub__(self,other):
		return Vector([self.data[0]-other.data[0],
					   self.data[1]-other.data[1],
					   self.data[2]-other.data[2]])
	def __getitem__(self,index):
		return self.data[index]
	def scale(self,scalar):
		return Vector([scalar*self.data[0],
					   scalar*self.data[1],
					   scalar*self.data[2]])
	def norm(self):
		return math.sqrt(self.data[0]*self.data[0]+
						 self.data[1]*self.data[1]+
						 self.data[2]*self.data[2])
	def hat(self):
		return self.scale(1.0/self.norm())
	def dot(self,other):
		return (self.data[0]*other.data[0]+
				self.data[1]*other.data[1]+
				self.data[2]*other.data[2])
	def cross(self,other):
		return Vector([self.data[1]*other.data[2]
					  -self.data[2]*other.data[1],
					   self.data[2]*other.data[0]
					  -self.data[0]*other.data[2],
					   self.data[0]*other.data[1]
					  -self.data[1]*other.data[0]])
	def get_x(self):
		return self.data[0]
	def get_y(self):
		return self.data[1]
	def get_z(self):
		return self.data[2]

# define crucial unit vectors
ehat= [Vector([1.0,0.0,0.0]),
	   Vector([0.0,1.0,0.0]),
	   Vector([0.0,0.0,1.0])]
Equinox= ehat[0]
CelestialPole= ehat[2]
EclipticAngle= 23.43928*piover180 # radians
EclipticPole= (ehat[2].scale(math.cos(EclipticAngle))
			  -ehat[1].scale(math.sin(EclipticAngle)))

# convert RA, Dec pair (in deg) to a three-vector
def RaDecToVector(aa,dd):
	arad= aa*piover180
	drad= dd*piover180
	return (
		ehat[0].scale(math.cos(drad)*math.cos(arad))
	   +ehat[1].scale(math.cos(drad)*math.sin(arad))
	   +ehat[2].scale(math.sin(drad)))

# returns (RA,Dec) in degrees
def VectorToRADec(vv):
	vv = vv.hat()
	dec = math.asin(vv.get_z())
	ra = math.atan2(vv.get_y(), vv.get_x())
	if ra < 0:
		ra += math.pi * 2
	return (ra*piunder180, dec*piunder180)
	

# from faint_motion.faint:
def mjdtodate(mjd):
	jd = mjdtojd(mjd)
	return jdtodate(jd)

def jdtodate(jd):
	unixtime = (jd - 2440587.5) * 86400. # in seconds
	return datetime.datetime.utcfromtimestamp(unixtime)

def mjdtojd(mjd):
	return mjd + 2400000.5

def timedeltatodays(dt):
	return dt.days + (dt.seconds + dt.microseconds/1e6)/86400.

def datetomjd(d):
	d0 = datetime.datetime(1858, 11, 17, 0, 0, 0)
	dt = d - d0
	# dt is a timedelta object.
	return timedeltatodays(dt)

# class for objects that are like planets, ie, objects with orbits around
#	  the Sun.
#	  name:		   planet name, of course
#	  elements:	   (a_AU, e, I_deg, L_deg, omega_deg, Omega_deg) at epoch
#	  elementsdot: time derivative of elements per century at epoch
#	  epoch:	   datetime UTC of epoch
class Planet:
	def __init__(self,name,elements,elementsdot,epoch,Vmag):
		self.Name= name
		self.ElementsAtEpoch= elements
		self.ElementDerivatives= elementsdot
		self.Epoch= epoch
		self.Vmag= Vmag
		self.ElementNames= ["a","e","I","L","omega","Omega"]
		self.ElementUnits= ["AU","","deg","deg","deg","deg"]
		self.Vector2= None
		self.DrawLabelEvenIfFaint= True
	def __repr__(self):
		return self.__str__()
	def __str__(self):
		return ('<Planet "%s": ' % self.Name
				+ ', '.join(['%s: %g %s' % (n,v,u)
							 for (n,v,u)
							 in zip(self.ElementNames, self.ElementsAtEpoch, self.ElementUnits)])
				+ ' at epoch ' + str(self.Epoch)
				+ '>')
	def copy(self):
		return Planet(self.Name, self.ElementsAtEpoch[:],
					  self.ElementDerivatives[:], self.Epoch,
					  self.Vmag)
	def set_elements(self, elems):
		self.ElementsAtEpoch = elems
	def set_element_derivs(self, delems):
		self.ElementDerivatives = delems
	def get_elements(self):
		return self.ElementsAtEpoch
	def get_element_derivs(self):
		return self.ElementDerivatives
	def Elements(self,UTC):
		deltat= UTC-self.Epoch
		deltaday= deltat.days+deltat.seconds/86400.0
		deltacentury= deltaday/36525.0
		elements= []
		ii= 0
		for element in self.ElementsAtEpoch:
			thiselement= (self.ElementsAtEpoch[ii]
						  +self.ElementDerivatives[ii]*deltacentury)
			if (self.ElementUnits[ii] == "deg"):
				while (thiselement > 180.0): thiselement-= 360.0
				while (thiselement < -180.0): thiselement+= 360.0
			elements.append(thiselement)
			ii+= 1
		return elements
	def RaDecFrom(self, viewpoint, UTC):
		vv = viewpoint.Vector3(UTC)
		vme = self.Vector3(UTC)
		return VectorToRADec(vme - vv)

	def get_PQmae(self, UTC):
		elements= self.Elements(UTC)
		semimajoraxis= elements[0]
		eccentricity= elements[1]
		inclination= elements[2]*piover180 # radians
		meananomaly= (elements[3]-elements[4])*piover180 # radians
		while (meananomaly > math.pi): meananomaly-= 2.0*math.pi
		while (meananomaly <= -math.pi): meananomaly+= 2.0*math.pi
		argumentperihelion= (elements[4]-elements[5])*piover180 # rad
		longascendingnode= elements[5]*piover180 # radians
		tmpydir= EclipticPole.cross(Equinox)
		ascendingnode= (Equinox.scale(math.cos(longascendingnode))
						+tmpydir.scale(math.sin(longascendingnode)))
		tmpydir= EclipticPole.cross(ascendingnode)
		orbitpole= (EclipticPole.scale(math.cos(inclination))
					-tmpydir.scale(math.sin(inclination)))
		tmpydir= orbitpole.cross(ascendingnode)
		periheliondir= (ascendingnode.scale(math.cos(argumentperihelion))
						+tmpydir.scale(math.sin(argumentperihelion)))
		tmpydir= orbitpole.cross(periheliondir)
		return (periheliondir, tmpydir, meananomaly*piunder180,
				semimajoraxis, eccentricity)

	def Vector3(self,UTC):
		vector3= Vector((0,0,0))
		if (self.ElementsAtEpoch[0] > 0):
			elements= self.Elements(UTC)
			semimajoraxis= elements[0]
			eccentricity= elements[1]
			inclination= elements[2]*piover180 # radians
			meananomaly= (elements[3]-elements[4])*piover180 # radians
			while (meananomaly > math.pi): meananomaly-= 2.0*math.pi
			while (meananomaly <= -math.pi): meananomaly+= 2.0*math.pi
			argumentperihelion= (elements[4]-elements[5])*piover180 # rad
			longascendingnode= elements[5]*piover180 # radians
			#print 'Long of asc =', rad2deg(longascendingnode)
			# y1
			tmpydir= EclipticPole.cross(Equinox)
			ascendingnode= (Equinox.scale(math.cos(longascendingnode))
							+tmpydir.scale(math.sin(longascendingnode)))
			#print 'Inclination =', rad2deg(inclination)
			# y2
			tmpydir= EclipticPole.cross(ascendingnode)
			orbitpole= (EclipticPole.scale(math.cos(inclination))
						-tmpydir.scale(math.sin(inclination)))
			#print 'Orbit pole =', orbitpole
			#print 'Arg of peri =', rad2deg(argumentperihelion)
			# y3
			tmpydir= orbitpole.cross(ascendingnode)
			periheliondir= (ascendingnode.scale(math.cos(argumentperihelion))
							+tmpydir.scale(math.sin(argumentperihelion)))
			#print 'Perihelion dir =', periheliondir
			#print 'Mean anomaly =', rad2deg(meananomaly), 'deg'
			# y4
			tmpydir= orbitpole.cross(periheliondir)
			#print 'Other perihelion dir =', tmpydir
			eccentricanomaly= meananomaly+eccentricity*math.sin(meananomaly)
			deltae= 100.0
			while (deltae > 1e-6):
				newmeananomaly= (eccentricanomaly
								 -eccentricity*math.sin(eccentricanomaly))
				deltam= meananomaly-newmeananomaly
				deltae= deltam/(1.0-eccentricity*math.cos(eccentricanomaly))
				eccentricanomaly+= deltae
			#print 'Eccentric anomaly =', rad2deg(eccentricanomaly)
			xx= semimajoraxis*(math.cos(eccentricanomaly)-eccentricity)
			yy= semimajoraxis*(math.sqrt(1-eccentricity*eccentricity)
							   *math.sin(eccentricanomaly))
			vector3= periheliondir.scale(xx)+tmpydir.scale(yy)
		return vector3
	def SetVector2(self,vector2):
		self.Vector2= vector2
		return True

# UTC for 2000 January 1.5
j2000= datetime.datetime(2000,1,1,12,0,0,0,tzinfo=None)

# make list of planets
def PlanetCatalog():
	planets= []
	planet= Planet(_('Sun'),(0,0,0,0,0,0),(0,0,0,0,0,0),j2000,-26.8)
	planets.append(planet)
	planet= Planet(_('Earth-Moon barycenter'),
				   (	1.00000261, 0.01671123,-0.00001531,
					   100.46457166,102.93768193,  0.0),
				   (  0.00000562,-0.00004392,-0.01294668,
					 35999.37244981, 0.32327364, 0.0),j2000,None)
	planets.append(planet)

	# MPC query results for "Other epoch = 2007-05-4.995":
	
	# Epoch (TT/JD) T (TT/JD)         Peri.     Node      Incl.       e             q             a      Epoch (TT)
	# 2454221.60833 2454224.9992012  24.25826 326.86764  19.11327 0.4324194     2.0531688      3.61740  2007/05/01.11

	e = 0.4324194
	a = 3.61740
	I = 19.11327
	omega = 24.25826
	Omega = 326.86764
	epoch = jdtodate(2454224.9992012)
	peridate = jdtodate(2454221.60833)
	periodyrs = sqrt(a**3)

	M = 0.
	L = M + Omega + omega
	pomega = Omega + omega

	planet= Planet(_('Comet Holmes'),
				   (a, e, I, L, pomega, Omega),
				   (0,0,0,
					360.*100./periodyrs, 0,0),
				   epoch, None)
	planets.append(planet)
	return planets



def old_plot_orbit(t0, emb, holmes):
	ex = []
	ey = []
	hx = []
	hy = []
	for i in range(0,350,10):
		t = t0 + datetime.timedelta(i)
		ve = emb.Vector3(t)
		ex.append(ve.get_x())
		ey.append(ve.get_y())
		vh = holmes.Vector3(t)
		hx.append(vh.get_x())
		hy.append(vh.get_y())
	clf()
	plot(array(ex), array(ey), 'r-',
		 array(hx), array(hy), 'b-',
		 array([ex[0]]), array([ey[0]]), 'ro',
		 array([hx[0]]), array([hy[0]]), 'bo')
	axhline(y=0)
	axvline(x=0)
	xlabel('x')
	ylabel('y')
	legend(('Earth', 'Holmes'))
	axis('equal')
	savefig('orbit.png')


def old_check_jpl():
	clf()
	ras1 = []
	decs1 = []
	ras2 = []
	decs2 = []
	times = []
	f = open('ephemerides.txt')
	for line in f:
		if line.startswith('#'):
			continue
		words = line.split()
		yr = int(words[0])
		mo = int(words[1])
		day = int(words[2])
		hr = int(words[3][:2])
		mn = int(words[3][2:4])
		sec = int(words[3][4:])
		t = datetime.datetime(yr, mo, day, hr, mn, sec)
		times.append(t)
		ra = float(words[4])
		ra *= 15.
		dec = float(words[5])
		ve = emb.Vector3(t)
		vh = holmes.Vector3(t)
		(ra2,dec2) = VectorToRADec(vh - ve)
		ras1.append(ra)
		decs1.append(dec)
		ras2.append(ra2)
		decs2.append(dec2)
	ras1 = array(ras1)
	ras2 = array(ras2)
	decs1 = array(decs1)
	decs2 = array(decs2)
	ras1 -= (ras1 > 180)*360
	ras2 -= (ras2 > 180)*360
	plot(ras1, decs1, 'ro-',
		 ras2, decs2, 'bx-')
	xlabel('RA (deg)')
	ylabel('Dec (deg)')
	legend(('JPL', 'me'))
	savefig('radec.png')


def plot_radec_zoom(imra, imdec, imrad, imt0,
					times, ras, decs):
	clf()
	I=find(imt0 > 0)
	medra = median(imra[I])
	meddec = median(imdec[I])
	medt = mjdtodate(median(imt0[I]))
	plot(ras, decs, 'r-',
		 [medra], [meddec], 'mo')
	for i in range(0, len(ras), 30):
		plot([ras[i]], [decs[i]], 'ro')
		text(ras[i],decs[i],str(times[i].date()))
	text(medra, meddec, str(medt.date()))
	legend(('comet', 'median'))
	xlabel('RA (deg)')
	ylabel('Dec (deg)')

	for (r,d,rad) in zip(imra,imdec,imrad):
		if rad > 5:
			continue
		c = Circle(xy=array([r,d]), radius=rad, alpha=0.2,
				   facecolor='b', edgecolor='k')
		c.set_clip_box(gca().bbox)
		gca().add_artist(c)

	nplotted = sum((imra > 40) * (imra < 80) *
				   (imdec > 30) * (imdec < 55) *
				   (imrad < 5))
	print 'N plotted:', nplotted

	axis((40,80,30,55))

	savefig('radec-zoom.png')


def plot_exif_times(images, ras, decs, times):
	ttrue = []
	terr = []
	ottrue = []
	oterr = []
	ot = []
	tt = []
	orad = []
	trad = []
	for (r,d,rad,t0,t1) in images:
		#if rad > 2:
		#	continue
		I = find(abs(ras - r)*cos(d) + abs(decs - d) < rad)
		## HACK - could match the crossover point.
		if len(I) == 0:
			continue
		mint = datetomjd(times[min(I)])
		maxt = datetomjd(times[max(I)])
		if t0:
			ottrue.append((mint + maxt)/2.)
			oterr.append((1 + maxt - mint)/2.)
			ot.append(datetomjd(t0))
			orad.append(rad)
		if t1:
			ttrue.append((mint + maxt)/2.)
			terr.append((1 + maxt - mint)/2.)
			tt.append(datetomjd(t1))
			trad.append(rad)
	ttrue = array(ttrue)
	terr = array(terr)
	ottrue = array(ottrue)
	oterr = array(oterr)
	ot = array(ot)
	tt = array(tt)
	orad = array(orad)
	trad = array(trad)

	clf()
	(tlo,thi) = (54398, 54440)
	#(tlo,thi) = (54300, 54600)
	plot([tlo,thi], [tlo,thi], '0.5')
	IT = find(trad < 2)
	IO = find(orad < 2)
	p2 = errorbar(ttrue[IT], tt[IT], xerr=terr[IT], fmt='o', c='r')
	p1 = errorbar(ottrue[IO], ot[IO], xerr=oterr[IO], fmt='o', c='b')
	legend((p1[0],p2[0]), ('Orig EXIF', 'EXIF'), loc='lower right')
	xlabel('True time (MJD)')
	ylabel('EXIF time (MJD)')
	xlim((tlo, thi))
	ylim((tlo, thi))
	savefig('times.png')

	clf()
	(tlo,thi) = (54398, 54440)
	plot([tlo,thi], [tlo,thi], '0.5')
	for (x,w,y) in zip(ttrue,terr,tt):
		r = Rectangle(xy=array([x,y]), width=w, height=1,
					  facecolor='r', alpha=1./w)
		r.set_clip_box(gca().bbox)
		gca().add_artist(r)
		r1 = r
	for (x,w,y) in zip(ottrue,oterr,ot):
		r = Rectangle(xy=array([x,y]), width=w, height=1,
					  facecolor='b', alpha=1./w)
		r.set_clip_box(gca().bbox)
		gca().add_artist(r)
		r2 = r
	legend((r2,r1), ('Orig EXIF', 'EXIF'), loc='lower right')
	xlabel('True time (MJD)')
	ylabel('EXIF time (MJD)')
	xlim((tlo, thi))
	ylim((tlo, thi))
	savefig('times2.png')


def plot_3space(emb, comet1, comet2, times, t, fn1, fn2):

	# true position of the comet at the timestamp times.
	vv = [comet1.Vector3(mjdtodate(tt)) for tt in t]
	truex = array([v.get_x() for v in vv])
	truey = array([v.get_y() for v in vv])
	truez = array([v.get_z() for v in vv])
	vv = [comet1.Vector3(tt) for tt in times]
	truexx = array([v.get_x() for v in vv])
	trueyy = array([v.get_y() for v in vv])
	truezz = array([v.get_z() for v in vv])

	vv = [comet2.Vector3(mjdtodate(tt)) for tt in t]
	myx = array([v.get_x() for v in vv])
	myy = array([v.get_y() for v in vv])
	myz = array([v.get_z() for v in vv])
	vv = [comet2.Vector3(tt) for tt in times]
	myxx = array([v.get_x() for v in vv])
	myyy = array([v.get_y() for v in vv])
	myzz = array([v.get_z() for v in vv])

	vv = [emb.Vector3(mjdtodate(tt)) for tt in t]
	earthx = array([v.get_x() for v in vv])
	earthy = array([v.get_y() for v in vv])
	earthz = array([v.get_z() for v in vv])
	vv = [emb.Vector3(tt) for tt in times]
	earthxx = array([v.get_x() for v in vv])
	earthyy = array([v.get_y() for v in vv])
	earthzz = array([v.get_z() for v in vv])

	#subplot(1, 2, 1)

	evendashes = [0,1,8,11]
	odddashes =  [0,11,8,1]

	h = plot(myx,   myz, 'go',
			 truex, truez, 'bo',
			 earthx, earthz, 'ro')
	plot(myxx,  myzz, 'g-',
		 truexx, truezz, 'b-',
		 earthxx, earthzz, 'r-')
	for (tx,tz,ex,ez,mx,mz) in zip(truex,truez,earthx,earthz,myx,myz):
		plot([ex,tx],[ez,tz], 'b--', dashes=evendashes)
		plot([ex,mx],[ez,mz], 'g--', dashes=odddashes)
		#plot([0,tx],[0,tz], 'b-')
		#plot([0,mx],[0,mz], 'g-')
	xlabel('x')
	ylabel('z')
	title('3-space positions')
	axis('equal')
	legend(h, ('Fit path', 'True path', 'Earth'), loc='lower left')
	savefig(fn1)

	clf()
	plot(myx,   myy, 'go',
		 truex, truey, 'bo',
		 earthx, earthy, 'ro')
	plot(myxx,  myyy, 'g-',
		 truexx, trueyy, 'b-',
		 earthxx, earthyy, 'r-')
	for (tx,ty,ex,ey,mx,my) in zip(truex,truey,earthx,earthy,myx,myy):
		plot([ex,tx],[ey,ty], 'b--', dashes=evendashes)
		plot([ex,mx],[ey,my], 'g--', dashes=odddashes)
		#plot([0,tx],[0,ty], 'b-')
		#plot([0,mx],[0,my], 'g-')
	xlabel('x')
	ylabel('y')
	title('3-space positions')
	axis('equal')
	legend(h, ('Fit path', 'True path', 'Earth'), loc='lower left')
	savefig(fn2)


def plot_fit(emb, myholmes, times, ras, decs,
			 imra, imdec, imrad, imcolors=None,
			 modelra=None, modeldec=None,
			 trajra=None, trajdec=None):
	if myholmes:
		myradecs = [myholmes.RaDecFrom(emb, tt) for tt in times]
		myras  = array([ra  for (ra,dec) in myradecs])
		mydecs = array([dec for (ra,dec) in myradecs])
	elif trajra and trajdec:
		myras = trajra
		mydecs = trajdec
	else:
		myradecs = []
		
	clf()
	xlabel('RA (deg)')
	ylabel('Dec (deg)')
	for i, (r,d,rad) in enumerate(zip(imra,imdec,imrad)):
		if rad > 5:
			continue
		if imcolors is not None:
			culur = imcolors[i]
		else:
			culur = 'g'
		c = Circle(xy=array([r,d]), radius=rad, alpha=0.2,
				   facecolor=culur, edgecolor='k')
		c.set_clip_box(gca().bbox)
		gca().add_artist(c)

		if modelra and modeldec:
			plot([r, modelra[i]], [d, modeldec[i]], 'k')

	#plot(imra, imdec, 'g.')
	h = plot(ras, decs, 'r-',
			 myras, mydecs, 'b-')
	plot(ras[::30], decs[::30], 'ro',
		 myras[::30], mydecs[::30], 'bo')
	legend((h[0], h[1]), ('True path', 'Fit'))
	axis((30,80,25,55))


# RA,Dec in degrees
def sample_images(n, ra, dec, maxradius):
	while True:
		J = numpy.random.permutation(len(ra))[:n]
		sdec = std(dec[J])
		print 'stddev(Dec) =', sdec
		if sdec > maxradius:
			continue
		sra = std(ra[J]) * math.cos(deg2rad(mean(dec[J])))
		print 'stddev(RA) =', sra
		if sra > maxradius:
			continue
		break
	return J


def initialize_orbit(ras, decs, times, emb):
	I = argsort(times)
	ras = ras[I]
	decs = decs[I]
	times = times[I]

	s = [(RaDecToVector(r,d) - emb.Vector3(mjdtodate(t))).hat()
		 for (r,d,t) in zip(ras, decs, times)]

	# Green uses 1-indexed arrays:
	sref = s[1]
	tref = times[1]

	R = emb.Vector3(mjdtodate(tref))
	## compute Rdot
	## FIXME - approximation...
	Rdot = (emb.Vector3(mjdtodate(tref + 1.)) - R)

	# mass of earth in kg / mass of sun in kg
	m = 5.9742e24 / 1.98892e30

	T3 = times[2] - tref
	T1 = tref - times[0]

	s1 = s[0]
	s3 = s[2]

	print 'T1', T1, 's1', s1
	print 'T3', T3, 's3', s3

	# eqn 7.38
	sdot = ( (sref - s1).scale(T3 / (T1 * (T1 + T3))) +
			 (s3 - sref).scale(T1 / (T3 * (T1 + T3))) )
	sdotdot = ( (s3 - sref).scale(2./(T3 * (T1 + T3))) -
				(sref - s1).scale(2./(T1 * (T1 + T3))) )

	## FIXME?
	k = 0.01720209895
	G = k**2
	m1 = 1.
	m2 = 0.

	# initial guess
	r = 2.
	for i in range(20):
		# eqn 7.43
		rho = (k**2 * ((1.+m)/R.norm()**3 - 1/r**3) *
			   sref.dot(sdot.cross(R)) /
			   sref.dot(sdot.cross(sdotdot)))
		# Triple product [a,b,c] = a . (b x c)
		# eqn 7.44
		r = sqrt(rho**2 + R.norm()**2 + 2.*rho*(R.dot(sref)))
		print 'rho=', rho, 'r=', r

	# eqn 7.45
	rhodot = ((k**2)/2. * (1./r**3 - (1.+m)/(R.norm()**3)) *
			  sref.dot(sdotdot.cross(R)) /
			  sref.dot(sdot.cross(sdotdot)))

	# eqn 7.40
	rdot = sref.scale(rhodot) + sdot.scale(rho) + Rdot
	# eqn 7.39
	rvec = sref.scale(rho) + R

	# in eqn 7.26, the vectors are in the ecliptic frame.
	y1 = EclipticPole.cross(Equinox)
	# rotate into the Ecliptic frame
	r1 = rvec.dot(Equinox)
	r2 = rvec.dot(y1)
	r3 = rvec.dot(EclipticPole)
	v1 = rdot.dot(Equinox)
	v2 = rdot.dot(y1)
	v3 = rdot.dot(EclipticPole)

	r = Vector([r1,r2,r3])
	v = Vector([v1,v2,v3])

	# eqn 6.7
	mu = G * (m1 + m2)

	# eqn 7.27
	a = mu * r.norm() / (2.*mu - r.norm() * v.norm()**2)

	# eqn 7.28; this is in ecliptic coords.
	h = r.cross(v)

	# eqn 7.30
	Omega = -atan2(h.get_y(), h.get_x())
	I = acos(h.get_x() / h.norm())

	# eqn 7.31
	e = sqrt(1. - h.norm()**2 / (mu*a))

	# eqn 7.32
	E = sign(r.dot(v)) * acos((a - r.norm()) / (a*e))

	# eqn 7.33
	M = E - e * sin(E)
	u = 2. * atan(sqrt((1+e)/(1-e)) * tan(E/2.))

	# eqn 7.35
	omega = (sign(r.get_z()) *
			 acos((r.get_x() * cos(Omega) + r.get_y() * sin(Omega)) / r.norm())) - u

	dL = 100.*360. / sqrt(a**3)

	return Planet('', 
				  (a, e, rad2deg(I), rad2deg(Omega + omega + M),
				   rad2deg(Omega + omega), rad2deg(Omega)),
				  (0, 0, 0, dL, 0, 0),
				  mjdtodate(tref), None)
	



### FIXME
# -this isn't actually a circular orbit, it's an ellipse:
#  I assumed the projection to (ecliptic) xy is a circle and the
#  z component is sinusoidal.  It's not too bad for moderate I.
def initialize_circular_orbit(ras, decs, times,
							  emb, cdist):
	cthetas = []
	czs = []
	#cxs = []
	#cys = []

	y1 = EclipticPole.cross(Equinox)
	for j in range(len(ras)):
		# the observed RA,Dec of the comet at time 't'
		ra = ras[j]
		dec = decs[j]
		t = times[j]
		tdate = mjdtodate(t)
		# the earth's position
		e = emb.Vector3(tdate)
		# the direction to the comet.
		v = RaDecToVector(ra, dec)
		# make the distance to the comet from sun be "cdist"
		ee = e.dot(e)
		vv = v.dot(v)
		ve = e.dot(v)
		L = (-2*ve + sqrt(4*ve**2 - 4*vv*(ee-cdist**2))) / (2.*vv)
		# comet position
		cpos = e + v.scale(L)
		# rotate into the Ecliptic frame
		ecx = cpos.dot(Equinox)
		ecy = cpos.dot(y1)
		ecz = cpos.dot(EclipticPole)
		# express in cylindrical coords
		#cxs.append(ecx)
		#cys.append(ecy)
		czs.append(ecz)
		cthetas.append(rad2deg(math.atan2(ecy, ecx)))

	ctheta = array(cthetas)
	cz = czs = array(czs)
	#cxs = array(cxs)
	#cys = array(cys)
	t = times
	theta = mean(ctheta)

	# orbital rate in degrees / day
	dthetadt = mean((ctheta - theta) / (t - mean(t)))
	print 'dtheta/dt =', dthetadt, '(target', 360./(365.25 * 6.88),')'
	# -> degrees / century
	dL = dthetadt * 365.25 * 100.

	z = mean(cz)
	dzdtheta = mean((cz - z) / deg2rad(ctheta - theta))
	thetaminusOmega = rad2deg(math.atan2(z, dzdtheta))
	Omega = theta - thetaminusOmega
	if Omega < 0:
		Omega += 360.
	z0 = z / math.sin(deg2rad(thetaminusOmega))
	I = rad2deg(math.asin(z0 / cdist))
	epoch = mjdtodate(mean(t))
	print 'I =', I, 'degrees'
	print 'z0 =', z0
	print 'Omega', Omega
	e5 = Omega
	e4 = e5
	e3 = (theta-Omega) + e4
	return Planet('', 
				  (cdist, 0, I, e3, e4, e5),
				  (0, 0, 0, dL, 0, 0),
				  epoch, None)


pc = PlanetCatalog()

emb = pc[1]
holmes = pc[2]

print 'Holmes:', str(holmes)

#print 'At epoch:'
#holmes.Vector3(holmes.Epoch)
#print 'P=(%g,%g,%g)' % (0.976, -0.212, 0.055)
#print 'Q=(%g,%g,%g)' % (0.127,  0.749, 0.650)

images = pickle.load(open('images.pickle'))
imra  = array([im[0] for im in images])
imdec = array([im[1] for im in images])
imrad = array([im[2] for im in images])
imt0  = array([im[3] and datetomjd(im[3]) or 0 for im in images])
imt1  = array([im[4] and datetomjd(im[4]) or 0 for im in images])

times = [mjdtodate(54282 + i) for i in range(300)]
radecs = [holmes.RaDecFrom(emb, t) for t in times]
ras  = array([ra  for (ra,dec) in radecs])
decs = array([dec for (ra,dec) in radecs])

#plot_radec_zoom(imra, imdec, imrad, imt0, times, ras, decs)
#plot_exif_times(images, ras, decs, times)

I = find((imt0 > 0) * (imrad < 2))
allra  = imra[I]
alldec = imdec[I]
allt   = imt0[I]


### TEST
J = sample_images(3, allra, alldec, 1)
myholmes = initialize_orbit(allra[J], alldec[J], allt[J], emb)
myholmes.Name = 'my Holmes'
print 'initialize_orbit:', myholmes

plot_fit(emb, myholmes, times, ras, decs, [], [], [])
plot(allra[J], alldec[J], 'mo')
savefig('init1.png')
sys.exit(0)

# MAGIC (second) 5: 5 deg std for initialization points
J = sample_images(5, allra, alldec, 5)
# MAGIC 2: distance in AU to comet
myholmes = initialize_circular_orbit(allra[J], alldec[J], allt[J],
									 emb, 2.)

myholmes = initialize_circular_orbit(allra[J], alldec[J], allt[J],
									 emb, 2.)

myholmes.Name = 'my Holmes'

plot_fit(emb, myholmes, times, ras, decs, [], [], [])
plot(allra[J], alldec[J], 'mo')
savefig('init.png')

# Optimize it...

orig_myholmes = myholmes.copy()

import levmar

initra = mean(allra[J])
initdec = mean(alldec[J])

cut1 = (imt0 > 0)
cut2 = ((imt0 == 0) * (imt1 > 0))
cut3 = (((imra - initra)*math.cos(deg2rad(initdec)))**2 + (imdec - initdec)**2 < 20**2)
imdates = imt0[:]
J2 = find(cut2)
imdates[J2] = imt1[J2]
J = find(cut1 * cut3)

steps = 0
costmargin = 0.01

def rolloff(x):
	if abs(x) < 1:
		return x
	elif abs(x) < 2:
		return sign(x)*sqrt(abs(x))
	else:
		return sign(x)*sqrt(2.)


def PQmaetovector(P,Q,m0,a,e,epoch,t):
	# degrees/day
	dmdt = 360. / (sqrt(a**3) * 365.25)
	# days
	dt = (t - epoch)
	# degrees
	m = m0 + dmdt * dt
	# -> radians
	m *= piover180
	meananomaly = m
	eccentricity = e
	# copy-paste from Planet.Vector3()
	while (meananomaly > math.pi): meananomaly -= 2.0*math.pi
	while (meananomaly <= -math.pi): meananomaly += 2.0*math.pi
	# ...
	eccentricanomaly = meananomaly+eccentricity*math.sin(meananomaly)
	deltae= 100.0
	while (deltae > 1e-6):
		newmeananomaly = (eccentricanomaly
						  -eccentricity*math.sin(eccentricanomaly))
		deltam = meananomaly-newmeananomaly
		deltae = deltam/(1.0-eccentricity*math.cos(eccentricanomaly))
		eccentricanomaly += deltae
	xx = a*(math.cos(eccentricanomaly)-eccentricity)
	yy = a*(math.sqrt(1-eccentricity*eccentricity)
			*math.sin(eccentricanomaly))
	vc = P.scale(xx) + Q.scale(yy)
	return vc

def PQmaetoRADec(P,Q,m0,a,e,epoch,t,viewpoint):
	vc = PQmaetovector(P,Q,m0,a,e,epoch,t)
	vv = viewpoint.Vector3(mjdtodate(t))
	return VectorToRADec(vc - vv)

bestcost = 1e100

def PQobjective(params, target, extra):
	global steps
	
	(p1,p2,p3, q1,q2,q3, m0, a, e) = params
	(epoch,) = extra
	P = Vector((p1,p2,p3))
	P = P.hat()
	Q = Vector((q1,q2,q3))
	Q = Q.hat()

	costs = []
	mra = []
	mdec = []
	for j in J:
		(modelra, modeldec) = PQmaetoRADec(P,Q,m0,a,e,epoch,
										   imdates[j], emb)
		mra.append(modelra)
		mdec.append(modeldec)
		dra = (imra[j] - modelra)*math.cos(deg2rad(imdec[j])) / imrad[j]
		ddec = (imdec[j] - modeldec) / imrad[j]
		costs.append(rolloff(dra))
		costs.append(rolloff(ddec))
	cost = sqrt(sum([c**2 for c in costs]))
	print 'cost:', cost

	global bestcost
	if cost+costmargin < bestcost:
		culurs = []
		for c in costs:
			if c**2 > 2:
				culurs.append('r')
			else:
				culurs.append('g')

		tradec = []
		for tt in times:
			tradec.append(PQmaetoRADec(P,Q,m0,a,e,epoch,datetomjd(tt),emb))
		tra = [ra for (ra,dec) in tradec]
		tdec = [dec for (ra,dec) in tradec]

		plot_fit(emb, None, times, ras, decs,
				 imra[J], imdec[J], imrad[J],
				 culurs, mra, mdec, tra, tdec)
		#title('cost=%.1f: ' % cost + 'a=%.2f, e=%.2f, I=%.1f, L=%.1f, om=%.1f, Om=%.1f' % params)
		fn = 'fit-%04i.png' % steps
		savefig(fn)
		print '                 wrote', fn
		steps += 1
		bestcost = cost

	return costs





from levmar_wrapper import *

epoch = myholmes.Epoch
allresults = []
allcosts = []

#for ang in range(0, 360, 30):
for ang in [0]:

	print 'Angle:', ang

	## What does the optimizer do at the true trajectory parameters?
	epoch = holmes.Epoch
	(P,Q,m,a,e) = holmes.get_PQmae(epoch)

	#(P,Q,m,a,e) = myholmes.get_PQmae(epoch)
	#e = 0.05
	pp = (P.scale(math.cos(deg2rad(ang))) +
		  Q.scale(math.sin(deg2rad(ang))))
	qq = (P.scale(-math.sin(deg2rad(ang))) +
		  Q.scale(math.cos(deg2rad(ang))))
	P = pp
	Q = qq
	m -= ang
	(p1,p2,p3) = (P.get_x(), P.get_y(), P.get_z())
	(q1,q2,q3) = (Q.get_x(), Q.get_y(), Q.get_z())

	params0 = (p1,p2,p3,q1,q2,q3,m,a,e)
	print 'True params:', params0
	print 'Initial cost:', sqrt(sum([c**2 for c in PQobjective(params0, None, (datetomjd(epoch),))]))
	print 'Best cost:', 2.3317174971

	# fit-* plots
	bestcost = 1e100

	da = 1.
	de = 0.5
	dm = 60.
	dp = 0.5

	opt = Optimizer(PQobjective,
					measurements=(0,) * 2*len(J),
					extra=(datetomjd(epoch),),
					scale=(dp,dp,dp, dp,dp,dp, dm, da, de))
	opt.jacobian_delta = 1e-6

	for optstep in range(1):
	#for optstep in range(5):

		opt.set_low((p1-dp,p2-dp,p3-dp,
					 q1-dp,q2-dp,q3-dp,
					 m-dm, max(0, a-da), max(0, e-de)))
		opt.set_high((p1+dp,p2+dp,p3+dp,
					  q1+dp,q2+dp,q3+dp,
					  m+dm, a+da, min(0.9, e+de)))
		opt.set_offset((p1,p2,p3,q1,q2,q3,m,a,e))

		result = opt.optimize()

		if not result:
			print 'WTF!'
			sys.exit(0)
		# renormalize P,Q.
		(p1,p2,p3, q1,q2,q3, m, a, e) = result
		P = Vector((p1,p2,p3))
		P = P.hat()
		Q = Vector((q1,q2,q3))
		Q = Q.hat()
		(p1,p2,p3) = (P.get_x(), P.get_y(), P.get_z())
		(q1,q2,q3) = (Q.get_x(), Q.get_y(), Q.get_z())
		result = (p1,p2,p3, q1,q2,q3, m, a, e)

		(trueP,trueQ,truem,truea,truee) = holmes.get_PQmae(epoch)
		print 'Result: P=(%.3f, %.3f, %.3f), Q=(%.3f, %.3f, %.3f), m=%.1f, a=%.2f, e=%.2f' % result
		print ('True  : P=(%.3f, %.3f, %.3f), Q=(%.3f, %.3f, %.3f), m=%.1f, a=%.2f, e=%.2f' %
			   (trueP.get_x(), trueP.get_y(), trueP.get_z(), trueQ.get_x(), trueQ.get_y(), trueQ.get_z(), truem, truea, truee))



		allresults.append(result)
		allcosts.append(bestcost)

		opt.set_initial(result)

		bestcost = 1e100
		PQobjective(result, [0.]*2*len(J), (datetomjd(epoch),))
		title('Angle %i, cost %.1f' % (ang, bestcost))
		savefig('fit-ang%03i-step%i.png' % (ang, optstep))

sys.exit(0)



def objective(params, target, comet):
	global steps

	(a,e,I,L,om,Om) = params
	comet.set_elements([a,e,I,L,om,Om])
	dL = 100.*360. / sqrt(a**3)
	comet.set_element_derivs([0,0,0,dL,0,0])

	costs = []
	mra = []
	mdec = []
	for j in J:
		t = mjdtodate(imdates[j])
		(modelra, modeldec) = comet.RaDecFrom(emb, t)
		mra.append(modelra)
		mdec.append(modeldec)
		dra = (imra[j] - modelra)*math.cos(deg2rad(imdec[j])) / imrad[j]
		ddec = (imdec[j] - modeldec) / imrad[j]
		costs.append(rolloff(dra))
		costs.append(rolloff(ddec))
	cost = sqrt(sum([c**2 for c in costs]))
	print 'cost:', cost

	global bestcost
	if cost+costmargin < bestcost:
		culurs = []
		for c in costs:
			if c**2 > 2:
				culurs.append('r')
			else:
				culurs.append('g')
				
		plot_fit(emb, comet, times, ras, decs,
				 imra[J], imdec[J], imrad[J],
				 culurs, mra, mdec)
		title('cost=%.1f: ' % cost + 'a=%.2f, e=%.2f, I=%.1f, L=%.1f, om=%.1f, Om=%.1f' % params)
		fn = 'fit-%04i.png' % steps
		savefig(fn)
		print '                 wrote', fn
		steps += 1
		bestcost = cost

	return costs


for optstep in range(10):

	(a,e,I,L,om,Om) = myholmes.get_elements()
	dL = 100.*360. / sqrt(a**3)

	dela = 1.
	dele = 0.5
	delI = 10.
	delL = 30.
	delom = 30.
	delOm = 30.

	opts = list(levmar.DEFAULT_OPTS) + [levmar.DIFF_DELTA]
	# ||J'e||_inf
	opts[1] = 1e-10
	# Dp, change in parameter values.
	opts[2] = 1e-6 # milliarcsec ~= 1e-6 deg, and RA,Dec are the best measured quantities
	# 1e-8 is more than enough precision, even for RA,Dec.
	# ||e||_2, norm of error.
	# opts[3] = 1e-3

	iterations = 10

	(result, iterations, run_info) = levmar.ddif_bc(
		objective, (a,e,I,L,om,Om), [0.]*2*len(J),
		[max(0, a-dela), max(0, e-dele),
		 I-delI, L-delL, om-delom, Om-delOm],
		[a+dela, min(0.95, e+dele),
		 I+delI, L+delL, om+delom, Om+delOm],
		iterations, opts=opts, data=myholmes)

	print 'Run info:', run_info


params = result

# make plot.
bestcost = 1e100

mycosts = array(objective(params, None, myholmes))
savefig('fit-final.png')

mycost = sqrt(sum([c**2 for c in mycosts]))
print 'My params:', params
print 'My cost:', mycost
print '(%i elements)' % len(mycosts)

trueparams = holmes.get_elements()
truecosts = array(objective(trueparams, None, holmes))
truecost = sqrt(sum([c**2 for c in truecosts]))
print 'True params:', trueparams
print 'True cost:', truecost
print '(%i elements)' % len(truecosts)

xxtimes = [mjdtodate(54000 + i) for i in range(3000)]
clf()
plot_3space(emb, holmes, myholmes, xxtimes, imdates[J],
			'3space-xz.png', '3space-xy.png')

clf()
bins = arange(0, 3.3, 0.1)
subplot(2,1,1)
hist(mycosts, bins)
xlabel('my costs')
subplot(2,1,2)
hist(truecosts, bins)
xlabel('true costs')
savefig('costs.png')

