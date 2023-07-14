from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    # other setup arguments here...
    include_dirs=[np.get_include()],
    ext_modules=cythonize("triplet_correlation_function.pyx"), name="ml_cy",
)
# After run: python setup.py build_ext --inplace
# Import like: import my_module