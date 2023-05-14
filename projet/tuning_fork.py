from ctypes import *

lib = "./python_functions.so"
functions = CDLL(lib)

functions.init_gmsh()
functions.get_k_freq_tuning_fork.restype = POINTER(c_double)
functions.get_k_freq_tuning_fork_n_layer.restype = POINTER(c_double)

def get_k_freq_tuning_fork(k, r1, r2, e, l, meshsize):
    res = functions.get_k_freq_tuning_fork(k, c_double(r1), c_double(r2), c_double(e), c_double(l), c_double(meshsize))
    ret = [res[i] for i in range(k)]
    functions.free_array(res)
    return ret

def get_k_freq_tuning_fork_n_layer(k, lh, wh, ws, btwhs, h, l, n, meshsize):
    ws_c = (c_double * len(ws)) (*ws)
    btwhs_c = (c_double * len(btwhs)) (*btwhs)
    h_c = (c_double * len(h)) (*h)
    l_c = (c_double * len(l)) (*l)
    res = functions.get_k_freq_tuning_fork_n_layer(k, c_double(lh), c_double(wh), ws_c, btwhs_c, h_c, l_c, n, c_double(meshsize))
    array = [res[i] for i in range(k)]
    functions.free_array(res)
    return array