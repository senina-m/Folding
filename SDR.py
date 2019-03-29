import scipy.optimize as sc
import numpy as np
import math

'''
Жалкая попытка повторить алгоритм SDR
попытка в разработке
'''


def fun(D, W, lambd):
    n = len(D[0])
    list = np.zeros(int(n * (n + 1) / 2))
    x = -1 / (n + math.sqrt(n))
    y = -1 / math.sqrt(n)
    V = np.concatenate(y * np.ones(1, n - 1), x * np.ones(n - 1) + np.eye(n - 1))
    e = np.ones(n, 1)

    G = np.zeros(n - 1, n - 1)
    B = V * G * V.T()
    E = np.diag(B) * e.T() + e * np.diag(B).T() - 2 * B

    def ineq_constraint(x):
        return x

    constraints = {'type': 'ineq', 'fun': ineq_constraint}

    def f(list):
        k = 0
        for i in range(n):
            for j in range(n):
                G[i][j] = list[k]
                G[j][i] = list[k]
            k = k + 1
        return np.trace(G) - lambd * np.linalg.norm(np.multiply(W, (E - D)))

    M = sc.minimize(f, list, method='nelder-mead', constraints=constraints, options={'xtol': 1e-8, 'disp': True})
    U, S, V = np.linalg.svd(B)
    EDM = np.diag(B) * e.T() + e * np.diag(B).T() - 2 * B
    X = math.sqrt(S) * V.T()
    return X
