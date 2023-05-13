import scipy
import math

from tuning_fork import get_k_freq_tuning_fork

def area_tuning_fork(r1, r2, e, l):
    return (r2 - r1) * e + 2 * (r2 - r1) * l + math.pi * (r2 ** 2 - r1 ** 2) / 4

def simple_mae(pred, target):
    return abs(pred - target)

def mae(preds, targets):
    return sum(abs(pred - target) for pred, target in zip(preds, targets))

def simple_loss_pred_g6(r1, r2, e, l):
    res = get_k_freq_tuning_fork(1, r1, r2, e, l)
    return simple_mae(res[0], 1567.98)

def minimize_simple_g6(r1, r2, e, l):
    res = scipy.optimize.minimize(lambda x: simple_loss_pred_g6(r1, *x), [r2, e, l], bounds=[(1e-3, 1.), (1e-3, 1.), (1e-3, 1.)], method="Nelder-Mead")
    return [r1] + list(res.x), res.fun

if __name__ == "__main__":
    r1, r2, e, l = 6e-3, 11e-3, 38e-3, 82e-3
    res, error = minimize_simple_g6(r1, r2, e, l)
    print(f"Dimension pour un simple diapason produisant la note G6, r1={res[0]}, r2={res[1]}, e={res[2]}, l={res[3]}, with an absolute error of {error}")