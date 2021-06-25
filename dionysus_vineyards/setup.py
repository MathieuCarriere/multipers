from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
extensions = [Extension('dionysus_vineyards', sources=['dionysus_vineyards.pyx'], language='c++')]
setup(name='dionysus_vineyards', ext_modules=cythonize(extensions), include_dirs=['.'])
