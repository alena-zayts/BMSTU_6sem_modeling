from math import exp
import numpy as np

a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1

a2 = 2.049
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1

k0 = 1.0
l = 10  # длина стержня
T0 = 300  # температура обдувающего воздуха
R = 0.5  # радиус стержня


Fmax = 50  # амплитуда импульса (определяет крутизну фронтов)
tmax = 60  # время (определяет длительность импульса)


# коэффициент теплопроводимости материала стержня
def lambda_T(n, m):
    return a1 * (b1 + c1 * pow(T[m][n],  m1))

def kappa_T(n, m):
    return lambda_T(n, m)

# коэффициент теплоемкости материала стержня
def c_T(n, m):
    return a2 + b2 * pow(T[m][n], m2) - c2 / (T[m][n] ** 2)


# коэффициент оптического поглощения материала стержня
def k_T(n, m):
    return k0 * pow(T[m][n] / 300, 2)

def exp_Tx(n, m):
    return exp(-k_T(n, m) * x_array[n])


alpha0 = 0.05
alphaN = 0.01
d_for_alpha = (alphaN * l) / (alphaN - alpha0)
c_for_alpha = - alpha0 * d_for_alpha
# коэффициент теплоотдачи при обдуве
def alpha_x(n):
    return c_for_alpha / (x_array[n] - d_for_alpha)

#!?
# импульс излучения

# Если в настоящей работе задать поток постоянным, т.е. F0(t) =const, то будет
# происходить формирование температурного поля от начальной температуры T0
# до некоторого установившегося (стационарного) распределения T(x,t)
# Это поле в дальнейшем с течением времени меняться не будет. Это полезный факт для тестирования программы.

# Если после разогрева стержня положить поток F(t)=0, то будет происходить
# остывание, пока температура не выровняется по всей длине и не станет равной T0
def F0_t(m):
    return Fmax * t_array[m] * exp(-(t_array[m] / tmax - 1)) / tmax

















EPS = 1e-6

xmin = 0
xmax = l
h = 1e-3

tmin = 0
tmax = tmax
tau = 1


x_array = list(np.arange(xmin, xmax + EPS, h))
t_array = list(np.arange(tmin, tmax + EPS, tau))

N = len(x_array)
M = len(t_array)


# Простая аппроксимация
def func_plus_half(func, n, m=None):
    if m:
        return (func(n, m) + func(n + 1, m)) / 2
    else:
        return (func(n) + func(n + 1)) / 2

def func_minus_half(func, n, m=None):
    if m:
        return (func(n, m) + func(n - 1, m)) / 2
    else:
        return (func(n) + func(n - 1)) / 2




T = []
T.append([T0] * N)


# поступает m, и считаю ддля шага m, то есть начать с 1
def count_boundary_conditions(m):
    prevT = T[m - 1]
    K0 = (h * func_plus_half(c_T, 0, m - 1) / 8 +
          h * c_T(0, m - 1) / 4 +
          tau * func_plus_half(kappa_T, 0, m-1) / h +
          h * tau * func_plus_half(alpha_x, 0) / (4 * R) +
          h * tau * alpha_x(0) / (2 * R) +
          alpha_x(0) * tau)

    M0 = (h * func_plus_half(c_T, 0, m - 1) / 8 -
          tau * func_plus_half(kappa_T, 0, m-1) / h +
          h * tau * func_plus_half(alpha_x, 0) / (4 * R))

    P0 = (h * func_plus_half(c_T, 0, m - 1) * (T[m - 1][0] + T[m - 1][1]) / 8 +
          h * c_T(0, m - 1) * T[m - 1][0] / 4 +
          alpha_x(0) * tau * T0 +
          h * tau * (
                F0_t(m - 1) * (k_T(0, m - 1) * exp_func(0, m - 1) + func_plus_half(k, 0, m - 1) * func_plus_half(exp_Tx, 0, m - 1)) +
                2 * T0 * (func_plus_half(alpha_x, 0) + alpha_x(n)) / R
          )/ 4)



    P0 = h / 8 * func_plus_1_2(prevT[0], t, c) * (T[m - 1][0] + T[m - 1][1]) + \
                        h / 4 * c(prevT[0]) * prevT[0] + \
                        F0 * t + t * h / 8 * (3 * f(0) + f(h))

    KN = h / 8 * func_minus_1_2(prevT[-1], t, c) + h / 4 * c(prevT[-1]) + \
                        func_minus_1_2(prevT[-1], t, k) * t / h + t * alphaN + \
                        t * h / 8 * p(l - h / 2) + t * h / 4 * p(l)

    MN = h / 8 * func_minus_1_2(prevT[-1], t, c) - \
                        func_minus_1_2(prevT[-1], t, k) * t / h + \
                        t * h * p(l - h / 2) / 8

    PN = h / 8 * func_minus_1_2(prevT[-1], t, c) * (prevT[-1] + prevT[-2]) + \
                        h / 4 * c(prevT[-1]) * prevT[-1] + t * alphaN * T0 + \
                        t * h / 4 * (f(l) + f(l - h / 2))









#######
def A(T):
    return func_minus_1_2(T, t, k) * t / h

def D(T):
    return func_plus_1_2(T, t, k) * t / h

def B(x, T):
    return A(T) + D(T) + c(T) * h + p(x) * h * t

def F(x, T):
    return f(x) * h * t + c(T) * T * h

def run_through(prevT):
    K0 = h / 8 * func_plus_1_2(prevT[0], t, c) + h / 4 * c(prevT[0]) + \
                        func_plus_1_2(prevT[0], t, k) * t / h + \
                        t * h / 8 * p(h / 2) + t * h / 4 * p(0)

    M0 = h / 8 * func_plus_1_2(prevT[0], t, c) - \
                        func_plus_1_2(prevT[0], t, k) * t / h + \
                        t * h * p(h / 2) / 8

    P0 = h / 8 * func_plus_1_2(prevT[0], t, c) * (prevT[0] + prevT[1]) + \
                        h / 4 * c(prevT[0]) * prevT[0] + \
                        F0 * t + t * h / 8 * (3 * f(0) + f(h))

    KN = h / 8 * func_minus_1_2(prevT[-1], t, c) + h / 4 * c(prevT[-1]) + \
                        func_minus_1_2(prevT[-1], t, k) * t / h + t * alphaN + \
                        t * h / 8 * p(l - h / 2) + t * h / 4 * p(l)

    MN = h / 8 * func_minus_1_2(prevT[-1], t, c) - \
                        func_minus_1_2(prevT[-1], t, k) * t / h + \
                        t * h * p(l - h / 2) / 8

    PN = h / 8 * func_minus_1_2(prevT[-1], t, c) * (prevT[-1] + prevT[-2]) + \
                        h / 4 * c(prevT[-1]) * prevT[-1] + t * alphaN * T0 + \
                        t * h / 4 * (f(l) + f(l - h / 2))

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

    # Обратный ход
    y = [0] * (n + 1)
    y[n] = (PN - MN * eta[n]) / (KN + MN * eps[n])

    for i in range(n - 1, -1, -1):
        y[i] = eps[i + 1] * y[i + 1] + eta[i + 1]

    return y

def check_eps(T, newT):
    for i, j in zip(T, newT):
        if fabs((i - j) / j) > 1e-2:
            return True
    return False

def check_iter(T, newT):
    max = fabs((T[0] - newT[0]) / newT[0])
    for i, j in zip(T, newT):
        d = fabs(i - j) / j
        if d > max:
            max = d
    return max < 1

def iter_method():
    result = []
    n = int(l / h)
    T = [T0] * (n + 1)
    newT = [0] * (n + 1)
    ti = 0

    result.append(T)

    while (True):
        buf = T
        while True:
            newT = run_through(buf)
            if check_iter(buf, newT):
                break
            buf = newT

        result.append(newT)
        ti += t
        if (check_eps(T, newT) == False):
            break

        T = newT

    return result, ti


class SimpleIterations:
    @staticmethod
    def stop_iteration(T_array_cur, T_array_prev):
        if len(T_array_prev) != N or len(T_array_cur) != N:
            raise ValueError('NNNN')
        differences = [abs((T_array_cur[n] - T_array_prev[n]) / T_array_cur[n]) for n in range(N)]
        max_difference = max(differences)
        if max_difference < EPS:
            return True
        else:
            return False

    @staticmethod
    # начиная с 1 до M-1
    def solve_step_m(m):
        T_beg = T[m - 1]


    # коэффициенты берутся с предыдущего шага
    # то есть первый раз m=1, а считать буду из 0 коэффициентов
    # m = 1,...
    # хз как правильно, пока

    @staticmethod
    def count_boundary_conditions(m):
        def count_K0():
            sum1 = h * func_plus_half(c, 0, )

    # поступает m, мы будем искать для m+1
    def solve(self, m):
        y_beg_array = T[m]
        return
