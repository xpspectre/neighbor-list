from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = [
    Extension(
        name='get_neighbors_c',
        sources=['get_neighbors_c.pyx'],
        extra_compile_args=['-std=c++11']
    )
]

setup(ext_modules=cythonize(ext))
