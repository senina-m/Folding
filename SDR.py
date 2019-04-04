import scipy.optimize as sc
import numpy as np
import math

'''
Жалкая попытка повторить алгоритм SDR
попытка в разработке
'''
D = [[0, 1, 4, 9, 16], [1, 0, 1, 4, 9], [4, 1, 0, 1, 4], [9, 4, 1, 0, 1], [16, 9, 4, 1, 0]]


def fun(D, W, lambd):
    n = len(D[0])
    list = np.zeros(int(n * (n + 1) / 2))
    x = -1 / (n + math.sqrt(n))
    y = -1 / math.sqrt(n)
    V = np.concatenate((np.dot(y, np.ones((1, n - 1))), np.dot(x, np.ones((n - 1))) + np.eye(n - 1)))
    e = np.ones((n, 1))

    G = np.zeros((n - 1, n - 1))
    B = np.dot(np.dot(V, G), V.transpose())
    m= np.dot(np.diag(B), (e.transpose()))
    E = np.dot(np.diag(B), (e.transpose())) + np.dot(e, (np.diag(B).T())) - np.multiply(2, B)

    def ineq_constraint(p):
        return p

    constraints = {'type': 'ineq', 'fun': ineq_constraint}

    def f(arr):
        k = 0
        for i in range(n):
            for j in range(n):
                G[i][j] = arr[k]
                G[j][i] = arr[k]
            k = k + 1
        return np.trace(G) - np.multiply(lambd, np.linalg.norm(np.multiply(W, (E - D))))

    M = sc.minimize(f, list, method='nelder-mead', constraints=constraints,
                    options={'xtol': 1e-8, 'disp': True})
    U, S, V = np.linalg.svd(B)
    EDM = np.dot(np.diag(B), (e.transpose())) + np.dot(e, (np.diag(B).transpose())) - np.multiply(2, B)
    X = np.multiply(math.sqrt(S), V.transpose())
    return X


print(fun(D, np.ones((5, 5)), 10))
