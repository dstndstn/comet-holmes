
#include "Python.h"

#include <stdio.h>
#include <assert.h>
#include <sys/param.h>

// numpy - this should be in site-packages/numpy/core/include
#include "arrayobject.h"

#define GM 2.95912e-04  // AU^3/d^2
#define c_au_per_yr 63239.6717 // google says

static const double Equinox[] = { 1, 0, 0 };
static const double Solstice[] = { 0, 0.91748213226576936, 0.39777699402184802 };
static const double EclipticPole[] = { 0, -0.39777699402184802, 0.91748213226576936};

//#define square(x) ((x)*(x))

static double square(double x) {
	return x*x;
}

// r = a cross b
static void cross(const double* a, const double*b, double* r) {
	r[0] = a[1]*b[2] - a[2]*b[1];
	r[1] = a[2]*b[0] - a[0]*b[2];
	r[2] = a[0]*b[1] - a[1]*b[0];
}

static void orbital_vectors_from_elements(double a, double e, double I,
										  double Omega, double pomega, double M,
										  double* xhat, double* yhat) {
	double ascendingnodevector[3];
	double tmpydir[3];
	double zhat[3];
	double sinI;
	double cosp,sinp;
	ascendingnodevector[0] = cos(Omega);
	ascendingnodevector[1] = sin(Omega);
	ascendingnodevector[2] = 0.0;
	//cross(khat, ascendingnodevector, tmpydir);
	tmpydir[0] = -ascendingnodevector[1];
	tmpydir[1] =  ascendingnodevector[0];
	tmpydir[2] = 0;
	// zhat= cos(i) * khat - sin(i) * tmpydir
	sinI = sin(I);
	zhat[0] = -sinI*tmpydir[0];
	zhat[1] = -sinI*tmpydir[1];
	zhat[2] = cos(I) -sinI*tmpydir[2];
	cross(zhat, ascendingnodevector, tmpydir);
	// xhat= cos(pomega) * ascendingnodevector + sin(pomega) * tmpydir
	cosp = cos(pomega);
	sinp = sin(pomega);
	xhat[0] = cosp * ascendingnodevector[0] + sinp * tmpydir[0];
	xhat[1] = cosp * ascendingnodevector[1] + sinp * tmpydir[1];
	xhat[2] = cosp * ascendingnodevector[2] + sinp * tmpydir[2];
	cross(zhat, xhat, yhat);
}

static double mean_anomaly_from_eccentric_anomaly(double E, double e) {
	return (E - e * sin(E));
}

static double eccentric_anomaly_from_mean_anomaly(double M, double e) {
	double E = M + e * sin(M);
	int iteration = 0;
	double deltaM = 100.0;

	int maximum_iteration = 100;
	double tolerance = 1e-15;

	while ((iteration < maximum_iteration) 
		   && (fabs(deltaM) > tolerance)) {
		deltaM = (M - mean_anomaly_from_eccentric_anomaly(E, e));
		E = E + deltaM / (1. - e * cos(E));
		iteration++;
	}
	return E;
}

static void observed_position_at_times(double a, double e, double I,
									   double Omega, double pomega, double M0,
									   double t0,
									   double* times, double* observer, int N, double* cometxyz) {
	double dMdt, b;
	int i;

	// aka 'meanfrequency'
	dMdt = sqrt(GM / (a*a*a));
	// semi-minor axis?
	b = a*sqrt(1. - e*e);

	for (i=0; i<N; i++) {
		int j,k;
		double r, dM=0, lastdM=0;
		double dx[3];
		double cometx[3];
		double M;
		double xhat[3], yhat[3];
		M = M0 + dMdt * (times[i] - t0);
		orbital_vectors_from_elements(a, e, I, Omega, pomega, M, xhat, yhat);
		// light-time correction loop; usually takes one step.
		for (j=0; j<100; j++) {
			double cosE, sinE, dEdt, c, d;
			double E = eccentric_anomaly_from_mean_anomaly(M-dM, e);
			cosE = cos(E);
			sinE = sin(E);
			dEdt = 1.0 / (1.0 - e * cosE) * dMdt;
			c = a * (cosE - e);
			d = b * sinE;
			r = 0;
			for (k=0; k<3; k++) {
				// FIXME -- could eliminate cometx[]
				cometx[k] = c * xhat[k] + d * yhat[k];
				dx[k] = cometx[k] - observer[3*i+k];
				r += dx[k]*dx[k];
			}
			// correct for light travel time delay.
			r = sqrt(r);
			dM = r / c_au_per_yr * dMdt;
			if (fabs(dM - lastdM) < 1e-12)
				break;
		}
		for (k=0; k<3; k++)
			dx[k] /= r;
		for (k=0; k<3; k++)
			cometxyz[3*i + k] = dx[0]*Equinox[k] + dx[1]*Solstice[k] + dx[2]*EclipticPole[k];
	}
}

static PyObject* comet_lnlikelihood(PyObject* self, PyObject* args,
									int do_chord) {
    int N, K;
    int i,j,k;
	double a,e,I,Omega,pomega,M0, t0;
    PyArrayObject* times;
    PyArrayObject* obs;
    PyArrayObject* imgxyz;
    PyArrayObject* imgrad;

    PyArrayObject* gQs;

	//PyArrayObject* tbests;

	double* cometxyz;
	double* obsxyz;
	double lnl;

    if (!PyArg_ParseTuple(args, "dddddddO!O!O!O!O!",
						  &a, &e, &I, &Omega, &pomega, &M0, &t0,
						  &PyArray_Type, &times,
						  &PyArray_Type, &obs,
						  &PyArray_Type, &imgxyz,
						  &PyArray_Type, &imgrad,
						  &PyArray_Type, &gQs)) {
		PyErr_SetString(PyExc_ValueError, "args: (a,e,i,Omega,pomega,M0, t0, times, observer_xyz, image_xyz, image_rad, gQ)");
		return NULL;
	}

	// FIXME -- check gQs.

    if (PyArray_NDIM(times) != 1) {
        PyErr_SetString(PyExc_ValueError, "time array must be one-dimensional");
        return NULL;
    }
    if (PyArray_NDIM(obs) != 2) {
        PyErr_SetString(PyExc_ValueError, "observer position array must be two-dimensional");
        return NULL;
    }
    if (PyArray_NDIM(imgxyz) != 2) {
        PyErr_SetString(PyExc_ValueError, "image position array must be two-dimensional");
        return NULL;
    }
    if (PyArray_NDIM(imgrad) != 1) {
        PyErr_SetString(PyExc_ValueError, "image radius array must be one-dimensional");
        return NULL;
    }
    if (PyArray_TYPE(times) != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "'times' array must contain doubles");
        return NULL;
    }
    if (PyArray_TYPE(obs) != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "'observer' array must contain doubles");
        return NULL;
    }
    if (PyArray_TYPE(imgxyz) != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "'imagexyz' array must contain doubles");
        return NULL;
    }
    if (PyArray_TYPE(imgrad) != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "'imgrad' array must contain doubles");
        return NULL;
    }
    if (PyArray_TYPE(gQs) != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError, "'gQ' array must contain doubles");
        return NULL;
    }

	//// FIXME -- what coord systems are 'observer' and 'imagexyz' in?
	// imagexyz should be from imageradec, so ecliptic unit vectors.
	// observer is in solar system coords?

    N = PyArray_DIM(times, 0);
	if (PyArray_DIM(obs, 0) != N) {
        PyErr_Format(PyExc_ValueError, "'obs' array must be Nx3 (N=%i)", N);
        return NULL;
	}
	if (PyArray_DIM(obs, 1) != 3) {
        PyErr_SetString(PyExc_ValueError, "'obs' array must be Nx3");
		return NULL;
	}

    K = PyArray_DIM(imgxyz, 0);
	if (PyArray_DIM(imgxyz, 1) != 3) {
        PyErr_Format(PyExc_ValueError, "'imgxyz' array must be Kx3 (K=%i)", K);
		return NULL;
	}
	if (PyArray_DIM(imgrad, 0) != K) {
        PyErr_Format(PyExc_ValueError, "'imgrad' array must be Kx1 (K=%i)", K);
		return NULL;
	}
	if (PyArray_DIM(gQs, 0) != K) {
        PyErr_Format(PyExc_ValueError, "'gQs' array must be Kx1 (K=%i)", K);
		return NULL;
	}

	/*
	 printf("a=%g, e=%g, i=%g, Omega=%g, omega=%g, M0=%g\n", a,e,I,Omega,pomega,M0);
	 printf("t0=%g\n", t0);
	 printf("N times: %i\n", N);
	 printf("K images: %i\n", K);
	 */
	/*{
		npy_intp nK = K;
		tbests = (PyArrayObject*)PyArray_SimpleNew(1, &nK, PyArray_DOUBLE);
	 }*/

	cometxyz = malloc(N * 3 * sizeof(double));

	// could avoid this if we knew the array was packed...
	// PyArray_GETCONTIGUOUS()...
	obsxyz = malloc(N * 3 * sizeof(double));
	for (i=0; i<N; i++) {
		//memcpy(obsxyz+i*3, PyArray_GETPTR2(obs, i, 0), 3*sizeof(double));
		for (j=0; j<3; j++)
			obsxyz[i*3+j] = *(double*)(PyArray_GETPTR2(obs, i, j));
	}

	observed_position_at_times(a, e, I, Omega, pomega, M0, t0,
							   PyArray_GETPTR1(times,0),
							   obsxyz, N, cometxyz);

	// evaluate max_i (Gaussian(imagexyz | cometxyz_i, imagerad))
	// or chord function...
	lnl = 0;
	for (k=0; k<K; k++) {
		double imxyz[3], imr;
		double minr2;
		double floorlnl;
		int imin;
		double tbest;
		double gQ;

		for (j=0; j<3; j++)
			imxyz[j] = *(double*)(PyArray_GETPTR2(imgxyz, k, j));
		imr = *(double*)(PyArray_GETPTR1(imgrad, k));

		gQ = *(double*)(PyArray_GETPTR1(gQs, k));

		minr2 = HUGE_VAL;
		for (i=0; i<N; i++) {
			double r2 = 0;
			for (j=0; j<3; j++)
				r2 += square(imxyz[j] - cometxyz[i*3+j]);
			if (r2 < minr2) {
				imin = i;
				minr2 = r2;
			}
		}
		tbest = *(double*)PyArray_GETPTR1(times, imin);

		assert(imin >= 1 && imin < (N-1));

		double r2a = 0;
		double r2b = 0;
		i = imin - 1;
		for (j=0; j<3; j++)
			r2a += square(imxyz[j] - cometxyz[i*3+j]);
		i = imin + 1;
		for (j=0; j<3; j++)
			r2b += square(imxyz[j] - cometxyz[i*3+j]);
		double *x0, *x1;
		double t0, t1;
		x0 = cometxyz + imin*3;
		t0 = tbest;
		if (r2a < r2b) {
			x1 = cometxyz + (imin-1)*3;
			t1 = *(double*)PyArray_GETPTR1(times, imin-1);
		} else {
			x1 = cometxyz + (imin+1)*3;
			t1 = *(double*)PyArray_GETPTR1(times, imin+1);
		}

		double dx[3];
		double rmin[3];
		double lendx2 = 0;
		double rmindotdx = 0;
		double frac;
		for (j=0; j<3; j++) {
			dx[j] = x1[j] - x0[j];
			rmin[j] = imxyz[j] - x0[j];
			lendx2 += dx[j]*dx[j];
			rmindotdx += rmin[j] * dx[j];
		}
		double velocity2 = lendx2 / square(t1 - t0);

		double rperp2;
		if (rmindotdx < 0) {
			// explode gently
			rperp2 = minr2;
		} else {
			rperp2 = 0;
			for (j=0; j<3; j++)
				rperp2 += square(rmin[j] - rmindotdx * dx[j] / lendx2);
		}

		if (do_chord) {
			double chord2 = 4 * MAX(0, imr*imr - rperp2);
			double timechord2 = chord2 / velocity2 ;
			lnl += log( 1.0 + gQ * sqrt(timechord2) );
		} else {
			double sigma = 0.5 * imr ;
			// totally made up by Hogg; ie, HACK //
			lnl += log( 1.0 + gQ * (imr / sqrt(velocity2))
						* exp(-0.5 * rperp2 / (sigma*sigma)) );
		}
	}
		
		/*
		 // DEBUG -- check difference between best and neighbours.
		 if (imin >= 1 && imin < (N-1)) {
		 double r2a = 0;
		 double r2b = 0;
		 i = imin - 1;
		 for (j=0; j<3; j++)
		 r2a += square(imxyz[j] - cometxyz[i*3+j]);
		 i = imin + 1;
		 for (j=0; j<3; j++)
		 r2b += square(imxyz[j] - cometxyz[i*3+j]);
		 double *x0, *x1;
		 double t0, t1;
		 x0 = cometxyz + imin*3;
		 t0 = tbest;
		 if (r2a < r2b) {
		 x1 = cometxyz + (imin-1)*3;
		 t1 = *(double*)PyArray_GETPTR1(times, imin-1);
		 } else {
		 x1 = cometxyz + (imin+1)*3;
		 t1 = *(double*)PyArray_GETPTR1(times, imin+1);
		 }
		 double dx[3];
		 double cdx[3];
		 double lendx2 = 0;
		 double cdxdotdx = 0;
		 double frac;
		 for (j=0; j<3; j++) {
		 dx[j] = x1[j] - x0[j];
		 cdx[j] = imxyz[j] - x0[j];
		 lendx2 += dx[j]*dx[j];
		 cdxdotdx += cdx[j] * dx[j];
		 }
		 // interpolate...
		 frac = MAX(0, cdxdotdx / lendx2);
		 minr2 -= frac * cdxdotdx;
		 tbest = t0 + (t1-t0) * frac;
		 //minr2 -= cdxdotdx*cdxdotdx / lendx2;
		 //tbest = t0 + (t1-t0) * cdxdotdx / lendx2;
		 //printf("interpolation fraction: %g\n", cdxdotdx / lendx2);

		 }

		 *(double*)PyArray_GETPTR1(tbests, k) = tbest;

		 // Don't let lnl go below 5-sigma (assume background model)
		 floorlnl = -(5*5 / 2.0);
		 lnl += MAX(floorlnl, -minr2 / (2.0*imr*imr));
		 */

	free(obsxyz);
	free(cometxyz);
    return Py_BuildValue("d", lnl);
    //return Py_BuildValue("(dO)", lnl, tbests);
}

static PyObject* comet_lnlikelihood_chord(PyObject* self, PyObject* args) {
	return comet_lnlikelihood(self, args, 1);
}

static PyObject* comet_lnlikelihood_gaussian(PyObject* self, PyObject* args) {
	return comet_lnlikelihood(self, args, 0);
}

static PyMethodDef cometLikelihoodMethods[] = {
    { "comet_lnlikelihood_chord", comet_lnlikelihood_chord, METH_VARARGS,
	  "stuff" },
    { "comet_lnlikelihood_gaussian", comet_lnlikelihood_gaussian, METH_VARARGS,
	  "stuff" },
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initcomet_likelihood(void) {
    Py_InitModule("comet_likelihood", cometLikelihoodMethods);
    import_array();
}


