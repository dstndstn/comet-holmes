import os.path
import datetime
import astrometry.util.defaults
import cPickle as pickle
from astrometry.util.starutil_numpy import *
from astrometry.util import EXIF
from astrometry.util.sip import *
from astrometry.util.file import *
from pylab import *
from numpy import *
from matplotlib.patches import Polygon

def wcs_polygon(wcs, refwcs, nsteps=1):
    radecs = wcs_polygon_radec(wcs, nsteps) 
    refxys = [refwcs.radec2pixelxy(*rd)
              for rd in radecs]
    return (array([rx for rx,ry in refxys]),
            array([ry for rx,ry in refxys]))

def wcs_polygon_radec(wcs, nsteps=1):
    xmin = 0.5
    xmax = wcs.imagew + 0.5
    x = concatenate((linspace(xmin, xmax, nsteps, endpoint=False),
                     array([xmax] * nsteps),
                     linspace(xmax, xmin, nsteps, endpoint=False),
                     array([xmin] * nsteps)))
    ymin = 0.5
    ymax = wcs.imageh + 0.5
    y = concatenate((array([ymin] * nsteps),
                     linspace(ymin, ymax, nsteps, endpoint=False),
                     array([ymax] * nsteps),
                     linspace(ymax, ymin, nsteps, endpoint=False)))
    #print 'x,y', x,y
    return [wcs.pixelxy2radec(xi,yi)
            for (xi,yi) in zip(x,y)]

def solid_angle(wcs):
    return arcsec2rad(arcsec2rad(wcs.imagew * wcs.imageh * (wcs.get_pixel_scale()**2)))

if __name__ == '__main__':
    dates = []
    ras = []
    decs = []
    radii = []
    imgfns = []
    wcsfns = []

    clf()
    
    refwcs = Tan()
    refwcs.crval[0] = 53
    refwcs.crval[1] = 50
    refwcs.cd[0] = 1.
    refwcs.cd[3] = 1.

    datafn = 'imagewcs.pickle'

    for i in range(1,1001):
        imgfn = 'images/yimage-%04i' % i
        wcsfn = 'images.wcs/yimage-%04i.wcs' % i
        if (not os.path.exists(imgfn)) or (not os.path.exists(wcsfn)):
            continue
        imgfns.append(imgfn)
        wcsfns.append(imgfn)

        print
        print imgfn
        exif = EXIF.process_file(open(imgfn), details=False)

        imgdatesrc = 'EXIF DateTimeOriginal'
        imgdate = exif.get(imgdatesrc)
        if not imgdate:
            imgdatesrc = 'Image DateTime'
            imgdate = exif.get(imgdatesrc)
        if not imgdate:
            imgdatesrc = None
            if not len(exif):
                print 'No exif'
            else:
                for k,v in exif.items():
                    print k,'=',v
        print imgdatesrc, imgdate

        if imgdate:
            date = datetime.datetime.strptime(str(imgdate), '%Y:%m:%d %H:%M:%S')
            dates.append(date)
        else:
            dates.append(None)

        wcs = Tan(filename=wcsfn)
        print 'Pixel scale:', wcs.get_pixel_scale(), 'arcsec/pixel'
        (r,d) = wcs.pixelxy2radec(wcs.imagew/2., wcs.imageh/2.)
        print 'RA,Dec center:', r,d
        ras.append(r)
        decs.append(d)

        radius = arcsec2deg(sqrt(wcs.imagew**2 + wcs.imageh**2)/2. * wcs.get_pixel_scale())
        radii.append(radius)

        omega = solid_angle(wcs)
        print 'omega', omega
        x,y = wcs_polygon(wcs, refwcs, nsteps=2)
        I = range(len(x)) + [0]
        
        p = Polygon(vstack((x[I], y[I])).T, alpha=max(0.00, 1e-5/omega),
                    edgecolor='none', facecolor='k')
        gca().add_artist(p)
        #plot(x[I], y[I], 'k-', lw=0.01/sqrt(omega))
        plot(x[I], y[I], 'k-', lw=0.2)

        if False:
            rds = wcs_polygon_radec(wcs, nsteps=100)
            I = range(len(rds)) + [0]
            plot([rds[i][0] for i in I],
                 sin(radians(array([rds[i][1] for i in I]))), 'r-')

    axis('scaled')
    ylim(-80,80)
    xlim(-80,80)
    savefig('rd0.png')

    axis('scaled')
    ylim(-20,20)
    xlim(-20,20)
    savefig('rd1.png')

    axis('scaled')
    ylim(-5,5)
    xlim(-5,5)
    savefig('rd2.png')

    write_file(pickle.dumps({'image-filenames': imgfns,
                             'wcs-filenames': wcsfns,
                             'dates': dates,
                             'center-ras': ras,
                             'center-decs': decs,
                             'radii': radii,
                             }), datafn)


    #I = argsort(dates)
    #med = I[len(dates)/2]
    #print 'median date:', dates[med]
    #print 'median ra,dec:', median(array(ras)), median(array(decs))
