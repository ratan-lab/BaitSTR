from distutils.core import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("fastq", ["fastq.pyx"]),
    Extension("select_STR_reads", ["select_STR_reads.pyx"]),
]

setup(version='1.0',  \
      ext_modules = cythonize(extensions),
     )
