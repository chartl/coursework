from distutils.core import setup, Extension

csmithwaterman = Extension('c_smithwaterman', ['./m260b/align/c_smithwaterman.cpp'], include_dirs=['/u/home/c/chartl/venv-base/lib/python2.7/site-packages/numpy/core/include/numpy'])
rlm = Extension('rlm', ['./m260b/assembly/rlm.cpp'], include_dirs=['/u/home/c/chartl/venv-base/lib/python2.7/site-packages/numpy/core/include/numpy'])

setup(name='m260b', version='0.0.1', description='m260b code', ext_modules=[csmithwaterman, rlm])
