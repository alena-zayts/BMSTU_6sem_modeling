import matplotlib.pyplot as plt
import numpy as np
from math import fabs, exp

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1
a2 = 2.049
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1
alpha0 = 0.05
alphaN = 0.01
l = 10
T0 = 300
R = 0.5
h = 0.01
t = 0.1
Fmax = 150
tmax = 10

period = 1
tu = 0.9

def k(T):
    return a1 * (b1 + c1 * T ** m1)
    '''K0 = 0.4
    KN = 0.1
    b = (KN * l) / (KN - K0)
    a = - K0 * b
    return a / (T - b)'''

def c(T):
    #return 0
    return a2 + b2 * T ** m2 - c2 / T ** 2

def alpha(x):
    d = (alphaN * l) / (alphaN - alpha0)
    c = - alpha0 * d
    return c / (x - d)

def p(x):
    return 2 * alpha(x) / R

def f(x):
    return 2 * alpha(x) * T0 / R

def func_plus_1_2(x, step, func):
    return (func(x) + func(x + step)) / 2

def func_minus_1_2(x, step, func):
    return (func(x) + func(x - step)) / 2

def F0(T):
    if (T % period <= tu):
        return Fmax / tmax * T * exp(-(T / tmax - 1))
    return 0
    #return Fmax / tmax * T * exp(-(T / tmax - 1))

def A(T):
    return func_minus_1_2(T, t, k) * t / h

def D(T):
    return func_plus_1_2(T, t, k) * t / h

def B(x, T):
    return A(T) + D(T) + c(T) * h + p(x) * h * t

def F(x, T):
    return f(x) * h * t + c(T) * T * h

def run_through(prevT, time):
    K0 = h / 8 * func_plus_1_2(prevT[0], t, c) + h / 4 * c(prevT[0]) + \
                        func_plus_1_2(prevT[0], t, k) * t / h + \
                        t * h / 8 * p(h / 2) + t * h / 4 * p(0)

    M0 = h / 8 * func_plus_1_2(prevT[0], t, c) - \
                        func_plus_1_2(prevT[0], t, k) * t / h + \
                        t * h * p(h / 2) / 8

    P0 = h / 8 * func_plus_1_2(prevT[0], t, c) * (prevT[0] + prevT[1]) + \
                        h / 4 * c(prevT[0]) * prevT[0] + \
                        F0(time) * t + t * h / 8 * (3 * f(0) + f(h))

    KN = h / 8 * func_minus_1_2(prevT[-1], t, c) + h / 4 * c(prevT[-1]) + \
                        func_minus_1_2(prevT[-1], t, k) * t / h + t * alphaN + \
                        t * h / 8 * p(l - h / 2) + t * h / 4 * p(l)

    MN = h / 8 * func_minus_1_2(prevT[-1], t, c) - \
                        func_minus_1_2(prevT[-1], t, k) * t / h + \
                        t * h * p(l - h / 2) / 8

    PN = h / 8 * func_minus_1_2(prevT[-1], t, c) * (prevT[-1] + prevT[-2]) + \
                        h / 4 * c(prevT[-1]) * prevT[-1] + t * alphaN * T0 + \
                        t * h / 4 * (f(l) + f(l - h / 2))

    #print(time, P0, F0(time))
    '''
    K0 = h / 8 * func_plus_1_2(prevT[0], t, c) + h / 4 * c(prevT[0]) + \
                        func_plus_1_2(0, h, k) * t / h + \
                        t * h / 8 * p(h / 2) + t * h / 4 * p(0)
    M0 = h / 8 * func_plus_1_2(prevT[0], t, c) - \
                        func_plus_1_2(0, h, k) * t / h + \
                        t * h * p(h / 2) / 8
    P0 = h / 8 * func_plus_1_2(prevT[0], t, c) * (prevT[0] + prevT[1]) + \
                        h / 4 * c(prevT[0]) * prevT[0] + \
                        F0(time) * t + t * h / 8 * (3 * f(0) + f(h))
    KN = h / 8 * func_minus_1_2(prevT[-1], t, c) + h / 4 * c(prevT[-1]) + \
                        func_minus_1_2(l, h, k) * t / h + t * alphaN + \
                        t * h / 8 * p(l - h / 2) + t * h / 4 * p(l)
    MN = h / 8 * func_minus_1_2(prevT[-1], t, c) - \
                        func_minus_1_2(l, h, k) * t / h + \
                        t * h * p(l - h / 2) / 8
    PN = h / 8 * func_minus_1_2(prevT[-1], t, c) * (prevT[-1] + prevT[-2]) + \
                        h / 4 * c(prevT[-1]) * prevT[-1] + t * alphaN * T0 + \
                        t * h / 4 * (f(l) + f(l - h / 2))
    '''
    # Прямой ход
    eps = [0, -M0 / K0]
    eta = [0, P0 / K0]

    x = h
    n = 1

    while (x + h < l):
        eps.append(D(prevT[n]) / (B(x, prevT[n]) - A(prevT[n]) * eps[n]))
        eta.append((F(x, prevT[n]) + A(prevT[n]) * eta[n]) / (B(x, prevT[n]) \
                                                        - A(prevT[n]) * eps[n]))
        n += 1
        x += h
    '''
    while (x + h < l):
        eps.append(D(x) / (B(x, x) - A(x) * eps[n]))
        eta.append((F(x, x) + A(x) * eta[n]) / (B(x, x) - A(x) * eps[n]))
        n += 1
        x += h
    '''

    # Обратный ход
    y = [0] * (n + 1)
    y[n] = (PN - MN * eta[n]) / (KN + MN * eps[n])

    for i in range(n - 1, -1, -1):
        y[i] = eps[i + 1] * y[i + 1] + eta[i + 1]

    return y

def check_eps(T, newT):
    for i, j in zip(T, newT):
        if fabs((i - j) / j) > 1e-4:
            return True
    return False

def check_iter(T, newT):
    max = fabs((T[0] - newT[0]) / newT[0])
    for i, j in zip(T, newT):
        d = fabs(i - j) / j
        if d > max:
            max = d
    return max < 1e-1

def iter_method():
    result = []
    n = int(l / h)
    T = [T0] * (n + 1)
    newT = [0] * (n + 1)
    ti = 0

    result.append(T)

    while (True):
        buf = T
        ti += t
        while True:
            newT = run_through(buf, ti)
            if check_iter(buf, newT):
                break
            buf = newT

        result.append(newT)
        if (check_eps(T, newT) == False):
            break

        T = newT

    return result, ti

def main():
    '''
    deltah = [1, 0.1, 0.01, 0.001]
    result = []
    global t
    for hi in deltah:
        t = hi
        res, ti = iter_method()
        x = [i for i in np.arange(0, l, h)]
        n = 0
        print(len(res))
        for temp in res:
            if (fabs(n-1) < 0.0001):
                print()
                print(t)
                print()
                result.append(temp[:-1])
            n += t
    print('    1    |   0.1   |   0.01  |  0.001  |')
    for i in range(len(result[0])):
        for j in range(len(result)):
            print(' %3.3f |' % result[j][i], end = '')
        print()
    deltah = [1, 0.1, 0.01, 0.001]
    result = []
    global h
    for hi in deltah:
        h = hi
        res, ti = iter_method()
        print(h)
        i = 0
        xfix = [temp[int(i / h)] for temp in res]
        result.append(xfix)
    print('    1    |   0.1   |   0.01  |  0.001  |')
    for i in range(len(result[0])):
        for j in range(len(result)):
            print(' %3.3f |' % result[j][i], end = '')
        print()
    arraya2 = [2.049, 5, 10, 15]
    arrayb2 = [0.000564, 0.001, 0.01, 0.1]
    result = []
    resulti = []
    global a2, b2
    for ai, bi in zip(arraya2, arrayb2):
        a2 = ai
        b2 = bi
        print(a2, b2)
        res, ti = iter_method()
        te = []
        i = 0
        while (i < ti):
            te.append(i)
            i += t
        i = 0
        xfix = [temp[int(i / h)] for temp in res]
        print(xfix[:-1])
        result.append(xfix[:-1])
        resulti.append(te)
    for res, teres in zip(result, resulti):
        plt.plot(teres, res)
    plt.xlabel("Время, c")
    plt.ylabel("Температура, K")
    plt.show()
    res, ti = iter_method()
    x = [i for i in np.arange(0, l, h)]
    n = 0
    for temp in res:
        if (n % period == 0):
            print(n)
            plt.plot(x, temp[:-1])
        n += 1
    plt.xlabel("Длина, см")
    plt.ylabel("Температура, K")
    plt.show()
    te = []
    i = 0
    while (i < ti):
        te.append(i)
        i += t
    i = 0
    xfix = [temp[int(i / h)] for temp in res]
    plt.plot(te, xfix[:-1])
    plt.xlabel("Время, c")
    plt.ylabel("Температура, K")
    plt.show()'''

    res, ti = iter_method()

    te = []
    i = 0
    while (i < ti):
        te.append(i)
        i += t
    i = 0
    xfix = [temp[int(i / h)] for temp in res]
    plt.plot(te, xfix[1:])
    plt.xlabel("Время, c")
    plt.ylabel("Температура, K")
    plt.show()

if __name__ == "__main__":
    main()
