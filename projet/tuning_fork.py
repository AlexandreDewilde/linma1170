from ctypes import *

lib = "./python_functions.so"
functions = CDLL(lib)

functions.init_gmsh()
functions.get_k_freq_tuning_fork.restype = POINTER(c_double)

def get_k_freq_tuning_fork(k, r1, r2, e, l):
    res = functions.get_k_freq_tuning_fork(k, c_double(r1), c_double(r2), c_double(e), c_double(l))
    return [res[i] for i in range(k)]

