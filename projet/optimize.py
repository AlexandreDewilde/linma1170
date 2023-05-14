import scipy
import math

from tuning_fork import get_k_freq_tuning_fork, get_k_freq_tuning_fork_n_layer

def area_tuning_fork(r1, r2, e, l):
    return (r2 - r1) * e + 2 * (r2 - r1) * l + math.pi * (r2 ** 2 - r1 ** 2) / 4

def simple_mae(pred, target):
    return abs(pred - target)

def mae(preds, targets):
    return sum(abs(pred - target) for pred, target in zip(preds, targets))

def simple_loss_pred_g6(r1, r2, e, l):
    res = get_k_freq_tuning_fork(1, r1, r2, e, l, 0.5)
    return simple_mae(res[0], 1567.98)

def minimize_simple_g6(r1, r2, e, l):
    res = scipy.optimize.minimize(lambda x: simple_loss_pred_g6(r1, *x), [r2, e, l], bounds=[(7e-3, 1e-1), (1e-3, 1e-1), (1e-3, 1e-1)], method="Nelder-mead")
    return [r1] + list(res.x), res.fun

def minimize_simple_g6_constant_area(area, r1, r2, e, l):
    res = scipy.optimize.minimize(lambda x: simple_loss_pred_g6(r1, *x) + abs(area_tuning_fork(r1, *x) - area)*10000000, [r2, e, l], bounds=[(7e-3, 1.), (1e-3, 30.), (1e-3, 1.)], method="Nelder-mead")
    return [r1] + list(res.x), res.fun

def constant_area(area, r1, r2, e, l):
    return area - area_tuning_fork(r1, r2, e, l)

def array_to_args(n, arr):
    ws = [arr[i] for i in range(2, n + 2)]
    btwhs = [arr[i] for i in range(n + 2, 2 * n + 2)]
    h = [arr[i] for i in range(2 * n + 2, 3 * n + 2)]
    l = [arr[i] for i in range(3 * n + 2, 4 * n + 2)]
    return arr[0], arr[1], ws, btwhs, h, l

def loss_pred_g6_n(n, x):
    # print(x)
    ret = get_k_freq_tuning_fork_n_layer(n, *array_to_args(n, x), n, 0.5)
    # print(ret)
    return mae(ret, [1567.98 * (i + 1) for i in range(n)])

def minimize_g6(lh, wh, ws, btwhs, h, l):
    classic = (1e-3, 1e-1)
    width = (1e-3, 1e-1)
    bounds = [classic, width, *[width]*len(ws), *[(0, 1e-2)]*len(btwhs), *[width]*len(h), *[classic]*len(h)]
    res = scipy.optimize.minimize(lambda x: loss_pred_g6_n(len(l), x), [lh] + [wh] + ws + btwhs + h + l, bounds=bounds)
    return list(res.x), res.fun

if __name__ == "__main__":
    r1, r2, e, l = 6e-3, 11e-3, 38e-3, 82e-3
    # res, error = minimize_simple_g6(r1, r2, e, l)
    # print(f"Dimension pour un simple diapason produisant la note G6, r1={res[0]}, r2={res[1]}, e={res[2]}, l={res[3]}, with an absolute error of {error}")
    
    # area = 1e-2
    res, error = minimize_simple_g6_constant_area(5e-4, r1, r2, e, l)
    print(area_tuning_fork(*res), error)
    print(get_k_freq_tuning_fork(1, *res, 0.5))
    # print(f"For a constant area {area}; the parameters are r1={res[0]}, r2={res[1]}, e={res[2]}, l={res[3]}, abs error is {error}")
    # res, error = minimize_g6(11e-3, 6e-3, [38e-3, 38e-3], [0,0], [32e-3, 32e-3], [82e-3, 82e-3])
    # print(res, error)