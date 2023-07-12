import os

__version__ = "0.0.1"

# set Python env variable to keep track of example data dir
DEBRaVE_dir = os.path.dirname(__file__)
DATADIR = os.path.join(DEBRaVE_dir, "TestSpectra/")

# Detect a valid CUDA environment
try:
    import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule

    cuda_ext = True
except:
    cuda_ext = False

try:
    from . import _kepler

    cext = True
except ImportError:
    cext = False
