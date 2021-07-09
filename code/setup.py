from distutils.core import setup, Extension
import numpy
import os.path

numpy_inc = os.path.join(numpy.get_include(), 'numpy')

c_module = Extension('comet_likelihood',
                     sources = ['comet_likelihood.c'],
                     include_dirs = [ numpy_inc, ],
                     extra_objects = [])

setup(name = 'comet likelihood',
      version = '1.0',
      description = '',
      author = 'Astrometry.net (Dustin Lang)',
      author_email = 'dstn@astro.princeton.edu',
      url = 'http://astrometry.net',
      py_modules = [ 'comet_likelihood' ],
      ext_modules = [c_module])

