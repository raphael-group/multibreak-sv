## setup.py code modified from 
## http://docs.cython.org/src/reference/compilation.html
## to distribute cython modules.
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "MBSV preprocessing app",
    ext_modules = cythonize('src/insDelHmm.pyx'),
    include_dirs=[numpy.get_include()],
)

