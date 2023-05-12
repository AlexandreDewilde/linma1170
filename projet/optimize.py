import subprocess
import scipy
import math

def f(r1, r2, e, l):
    # print(r1, r2, e, l)
    res = subprocess.getoutput(f"./execute 4 {r1} {r2} {e} {l} -v 0")
    print(res)
    return float(res)
    

r1,r2,e,l = 6e-3, 11e-3, 38e-3, 82e-3
r3,r4,l2 = 6e-3, 11e-3, 82e-3
d1 = r1
d2 = r2

def f2(e, l):
    res = subprocess.getoutput(f"./execute 4 {r1} {r2} {e} {l} -v 0")
    print(res)
    return float(res)

def f3(l2):
    res = subprocess.getoutput(f"./execute_2_layer {r1} {r2} {e} {l} {l2} -v 0")
    print(res)
    return float(res)

def f4(d1, d2, r1, r2, l, l2):
    res = subprocess.getoutput(f"./execute_2_layer {e} {d1} {d2} {r1} {r2} {l} {l2} -v 0")
    print(res)
    return float(res)

def area(r1, r2, e, l):
    return (r2 - r1) * e + 2 * (r2 - r1) * l + math.pi * (r2**2 - r1**2) / 4

target_area = area(r1,r2,e,l)
def con(x):
    return abs(area(r1, r2, x[0], x[1]) - target_area)

print(e, d1, d2, r1, r2, l, l2)
# print(r1, r2, e, l)
res = scipy.optimize.minimize(lambda x: f4(*x), [d1, d2, r1, r2, l, l2], bounds=[(1e-3, 1.), (1e-3, 1.), (5e-3, 1.), (5e-3, 1.), (1e-3, 1.), (1e-3, 1.)], options={"disp":True}, method="Nelder-mead")

# res = scipy.optimize.minimize(lambda x: f4(*x), [r2, e, l], bounds=[(7*1e-3, 1.), (1e-4, 1.), (1e-4, 1.)], options={"disp":True}, method="Nelder-Mead")
# d1, r1, l = res.x
print(res.x)
# print(f"{e} {d1} {d2} {r1} {r2} {l} {l2}")
