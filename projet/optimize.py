import subprocess
import scipy
import math

def f(r1, r2, e, l):
    # print(r1, r2, e, l)
    res = subprocess.getoutput(f"./execute 4 {r1} {r2} {e} {l} -v 0")
    print(res)
    return float(res)
    

r1,r2,e,l = 6e-3, 11e-3, 38e-3, 82e-3

def f2(e, l):
    print(r1, r2, e, l)
    res = subprocess.getoutput(f"./execute 4 {r1} {r2} {e} {l} -v 0")
    print(res)
    return float(res)

def area(r1, r2, e, l):
    return (r2 - r1) * e + 2 * (r2 - r1) * l + math.pi * (r2**2 - r1**2) / 4

target_area = area(r1,r2,e,l)
def con(x):
    return abs(area(r1, r2, x[0], x[1]) - target_area)

INF = 1.

m = 1e-2
print(r1, r2, e, l)
res = scipy.optimize.minimize(lambda x: f2(*x) + con(x), [e,l], bounds=[(1e-3, 1e-1), (1e-3, 1e-1)], options={"disp":True}, method="Nelder-Mead")
print(res)
