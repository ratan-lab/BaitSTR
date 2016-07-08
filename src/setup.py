from distutils.core import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension("fastq", ["fastq.pyx"]),
    Extension("select_STR_reads", ["select_STR_reads.pyx"]),
]

with open('VERSION',"r") as version_file:
    version = version_file.read().strip()

setup(version = version,  \
      ext_modules = cythonize(extensions),
     )
