# -*- coding: utf-8 -*-
import numpy as np


def pnleg(x, n):
    if n == 0:
        return 1
    else:
        p1, p2, p3 = 1, x, x
        for k in range(1, n):
            p3 = ((2 * k + 1) * x * p2 - k * p1) / (k + 1)
            p1, p2 = p2, p3
    return p3


def pnleg1(x, n):
    if n == 0:
        return 0, 1
    elif n == 1:
        return 1, x
    else:
        p1 = 0
        p2 = 0
        p11 = 0
        p21 = 0
        p3 = 0
        p31 = 0
        p1 = 1
        p2 = x
        p11 = 0
        p21 = 1

        for k in range(1, n):
            k2p1 = 2 * k + 1
            kp1 = 1 / (k + 1)
            p3 = (k2p1 * x * p2 - k * p1) * kp1
            p31 = (k2p1 * (x * p21 + p2) - k * p11) * kp1
            p11, p21 = p21, p31
            p1, p2 = p2, p3

        return p31, p3


def pnleg2(x, n):
    if n <= 2:
        return n - 1, n - 1, n - 1
    else:
        p1 = 1
        p2 = x
        p11 = 0
        p21 = 1
        p12 = 0
        p22 = 0
        for k in range(1, n):
            p3 = ((2 * k + 1) * x * p2 - k * p1) / (k + 1)
            p31 = ((2 * k + 1) * (x * p21 + p2) - k * p11) / (k + 1)
            p32 = ((2 * k + 1) * (x * p22 + p21 * 2) - k * p12) / (k + 1)
            p1, p2 = p2, p3
            p11, p21 = p21, p31
            p12, p22 = p22, p32
    return p32, p31, p3


def der2lgl(x, n):
    nm = n - 1
    d = np.zeros((n, n))
    for k in range(n):
        lnxl = pnleg(x[k], nm)
        if k != 0:
            d[0, k] = (
                (-1) ** nm / lnxl * (n * nm * (1 + x[k]) - 4) * 0.5 / (1 + x[k]) ** 2
            )
        for j in range(1, nm):
            if j != k:
                lnxj = pnleg(x[j], nm)
                d[j, k] = -2 * lnxj / (lnxl * (x[j] - x[k]) ** 2)
            else:
                d[j, j] = pnleg2(x[j], nm)[0] / (3 * lnxl)
        if k != nm:
            d[nm, k] = 1 / lnxl * (nm * n * (1 - x[k]) - 4) * 0.5 / (1 - x[k]) ** 2
    d[0, 0] = n * nm * (nm * nm + nm - 2) / 24
    d[nm, nm] = d[0, 0]
    return d


def derlgl(x, n):
    d = np.zeros((n, n))
    nm = n - 1
    for j in range(n):
        lnxj = pnleg(x[j], nm)
        for i in range(n):
            if i != j:
                lnxi = pnleg(x[i], nm)
                d[i, j] = lnxi / ((x[i] - x[j]) * lnxj)
    d[0, 0] = -0.25 * n * nm
    d[-1, -1] = -d[0, 0]
    return d


def jacobi_eval(x, n, alpha, beta):
    apb = alpha + beta
    ab2 = alpha * alpha - beta * beta
    if n == 0:
        return 1, 0
    else:
        p1, p2, pd1, pd2 = 1, 0, 0, 0
        p = (alpha - beta + (apb + 2) * x) * 0.5
        pd = 0.5 * (apb + 2) * np.ones_like(x)
        for k in range(1, n):
            k1 = k + 1
            k2 = k * 2
            k2ab = k2 + alpha + beta
            k2ab1 = k2ab + 1
            k2ab2 = k2ab1 + 1
            p2, p1 = p1, p
            pd2, pd1 = pd1, pd
            a1 = 2 * k1 * (k1 + apb) * k2ab
            a21 = k2ab1 * ab2
            a22 = k2ab2 * k2ab1 * k2ab
            a3 = 2 * (k + alpha) * (k + beta) * k2ab2
            p = ((a21 + a22 * x) * p1 - a3 * p2) / a1
            pd = (a22 * p1 + (a21 + a22 * x) * pd1 - a3 * pd2) / a1
    return p, pd


def jacobi_roots(n, alpha, beta):
    TOL = 1.0e-14
    KMAX = 15

    if n < 1:
        raise ValueError("n must be greater than 0")
    else:
        x = np.zeros(n)
        x0 = np.cos(np.pi / (2 * n))
        for j in range(n):
            diff = TOL + 1
            kiter = 0
            while kiter <= KMAX and diff >= TOL:
                p, pd = jacobi_eval(x0, n, alpha, beta)
                ss = np.sum(1 / (x0 - x[:j]))
                x1 = x0 - p / (pd - ss * p)
                diff = np.abs(x1 - x0)
                kiter += 1
                x0 = x1
            x0 = (x1 + np.cos((2 * (j + 2) - 1) * np.pi / (2 * n))) / 2.0
            x[j] = x1
    return np.sort(x)


def xwlgl(n):
    if n <= 1:
        raise ValueError("n must be greater than 1")
    elif n == 2:
        return np.array([-1, 1]), np.array([1, 1])
    else:
        x = np.zeros(n)
        w = np.zeros(n)
        nm = n - 1
        x[0] = -1
        x[-1] = 1
        x[1:-1] = jacobi_roots(nm - 1, 1, 1)
        pnlg = pnleg(x, nm)
        w = 2.0 / pnlg / pnlg / n / nm
    return x, w


def test_all_matlab():
    for i in range(2, 17):
        w_test = np.genfromtxt(f"test/wlgl{i:02d}.csv", delimiter=",")
        x_test = np.genfromtxt(f"test/xlgl{i:02d}.csv", delimiter=",")
        d_test = np.genfromtxt(f"test/dlgl{i:02d}.csv", delimiter=",")
        h_test = np.genfromtxt(f"test/hlgl{i:02d}.csv", delimiter=",")
        x, w = xwlgl(i)
        d = derlgl(x, i)
        h = der2lgl(x, i)

        print("Deg: ", i - 1)
        print(np.abs(w_test - w).max())
        print(np.abs(x_test - x).max())
        print(np.abs(d_test - d).max())
        print(np.abs(h_test - h).max())
        print("==================================")


def test_conv():
    for i in range(2, 17):
        x, w = xwlgl(i)
        d = derlgl(x, i)
        h = der2lgl(x, i)

        u = np.sin(x)
        u1 = d @ u
        u2 = h @ u

        print("Deg: ", i - 1)
        print(np.abs(u1 - np.cos(x)).max())
        print(np.abs(u2 + u).max())
        print("==================================")


def basis(x, n):
    xi, wi = xwlgl(n)
    coef = 1 / n / (n + 1)
    x_xi = np.subtract.outer(xi, x)
    d = np.empty_like(x_xi)
    for i in range(n):
        coef = np.prod(np.delete(xi, i) - xi[i])
        d[i] = np.prod(np.delete(x_xi, i, axis=0), axis=0) / coef
    return d


def basis_v2(x, n):
    xi, wi = xwlgl(n)
    V = np.power.outer(xi, np.arange(n))
    Vinv = np.linalg.inv(V)
    d = np.empty((n, x.size))
    for i in range(n):
        d[i] = np.polyval(Vinv[::-1, i], x)
    return d


def derlgl_v2(x, n):
    xi, wi = xwlgl(n)
    V = np.power.outer(xi, np.arange(n))
    Vinv = np.linalg.inv(V)
    d = np.empty((n, x.size))
    for i in range(n):
        d[i] = np.polyval(np.polyder(Vinv[::-1, i]), x)
    d2 = V[:, :-1] @ (Vinv[1:, :] * np.arange(1, n)[..., np.newaxis])
    return d.T


if __name__ == "__main__":
    # test_all_matlab()
    # test_conv()
    import time
    for n in range(3, 16):
        x, w = xwlgl(n)
        d = np.genfromtxt(f"test/dlgl{n:02d}.csv", delimiter=",")
    
        xx = np.linspace(-1, 1, 1000)
        b = basis(xx, n)
        b2 = basis_v2(xx, n)
    
        d2 = derlgl_v2(x, n)
    
        # import matplotlib.pyplot as plt
    
        # plt.plot(xx, b.T)
        # plt.vlines(x, -1, 1, color="k", linestyle="--")
        # plt.ylim([-1, 1])
    
        print(f"{np.max(np.abs(b2 - b)):.4e}")
        print(f"{np.max(np.abs(d2 - d)):.4e}")
