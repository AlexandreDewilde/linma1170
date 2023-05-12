import subprocess
import scipy
import math
from ctypes import *
lib = "./python_functions.so"
functions = CDLL(lib)


r1,r2,e,l = 6e-3, 11e-3, 38e-3, 82e-3
r3,r4,l2 = 6e-3, 11e-3, 82e-3
d1 = r1
d2 = r2

functions.get_k_freq_tuning_fork.restype = POINTER(c_double)
res = ((functions.get_k_freq_tuning_fork(2, c_double(r1), c_double(r2), c_double(e), c_double(l))))
print(res[0], res[1])

def area(r1, r2, e, l):
    return (r2 - r1) * e + 2 * (r2 - r1) * l + math.pi * (r2**2 - r1**2) / 4

target_area = area(r1,r2,e,l)
def con(x):
    return abs(area(r1, r2, x[0], x[1]) - target_area)

# print(e, d1, d2, r1, r2, l, l2)
# print(r1, r2, e, l)
# res = scipy.optimize.minimize(lambda x: f4(*x), [d1, d2, r1, r2, l, l2], bounds=[(1e-3, 1.), (1e-3, 1.), (5e-3, 1.), (5e-3, 1.), (1e-3, 1.), (1e-3, 1.)], options={"disp":True}, method="Nelder-mead")

# res = scipy.optimize.minimize(lambda x: f4(*x), [r2, e, l], bounds=[(7*1e-3, 1.), (1e-4, 1.), (1e-4, 1.)], options={"disp":True}, method="Nelder-Mead")
# d1, r1, l = res.x
# print(res.x)
# print(f"{e} {d1} {d2} {r1} {r2} {l} {l2}")
