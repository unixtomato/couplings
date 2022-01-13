import numpy as np
from scipy.optimize import root, fsolve

import timeit
import sys

# constants
l = 32
n_reno = 4
n_coup = 7

method = 'maj'
mults = np.array([2, 2, 2, 2, 2, 4, 4])


indices = np.array([(l // 2**r)**2 for r in range(n_reno)])

# load data
data = np.fromfile("plots/callen_l{}_{}.dat".format(l, method), dtype=np.int32)
data = data.reshape(-1, np.sum(indices) * n_coup)
total_size = data.shape[0]

data = np.split(data, np.cumsum(indices)[:-1] * n_coup,  axis=1)
data = [data[r].reshape(-1, indices[r], n_coup) for r in range(n_reno)]

# load corr from first expression
ops_data = np.fromfile("plots/ops_l{}_{}.dat".format(l, method), dtype=np.int32)
ops_data = ops_data.reshape(-1, n_reno, n_coup)
ops = np.mean(ops_data, axis=0) / indices[..., np.newaxis]

n_bins = 10
b_size = total_size // n_bins


def fun(coups, r, b):

    d = data[r][b*b_size:(b+1)*b_size, :, :]
    dd = np.dot(d, coups)
    f = np.mean(np.sum(np.tanh(dd[..., np.newaxis]) * d, axis=1), axis=0) / mults / indices[r]

    return f - ops[r]

def test():

    start_time = timeit.default_timer()
    
    # read in operators (correlation function)

    coups = np.array([0.44068679350977147, 0., 0., 0., 0., 0., 0.])

    for r in range(n_reno):
        bins = []
        for b in range(n_bins):
            sol = root(fun, coups, jac=False, args=(r,b))
            bins.append(sol.x)

            print(" %d/%d %d/%d" % (r, n_reno, b, n_bins),
                        (timeit.default_timer() - start_time), file=sys.stderr, end='\r')

        bins = np.array(bins).transpose()
        for c in range(n_coup):
            coups[c] = np.mean(bins[c]) # update coups to be used in the next renorm as init
            print(coups[c], np.std(bins[c])/np.sqrt(n_bins-1))


    end_time = timeit.default_timer()
    print("took", (end_time - start_time), "seconds")


if __name__ == '__main__':
    test()

