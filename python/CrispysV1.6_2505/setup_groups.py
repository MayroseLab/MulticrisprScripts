from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize('Covers.pyx'),  # accepts a glob pattern
)
