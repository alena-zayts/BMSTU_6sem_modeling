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
def lambda_T(T):
    return a1 * (b1 + c1 * pow(T,  m1))


# коэффициент теплоемкости материала стержня
def c_T(T):
    return a2 + b2 * pow(T, m2) - c2 / (T ** 2)


# коэффициент оптического поглощения материала стержня
def k_T(T):
    return k0 * pow(T / 300, 2)



alpha0 = 0.05
alphaN = 0.01
d_for_alpha = (alphaN * l) / (alphaN - alpha0)
c_for_alpha = - alpha0 * d_for_alpha
# коэффициент теплоотдачи при обдуве
def alpha_x(x):
    return c_for_alpha / (x - d_for_alpha)

#!?
# импульс излучения

# Если в настоящей работе задать поток постоянным, т.е. F0(t) =const, то будет
# происходить формирование температурного поля от начальной температуры T0
# до некоторого установившегося (стационарного) распределения T(x,t)
# Это поле в дальнейшем с течением времени меняться не будет. Это полезный факт для тестирования программы.

# Если после разогрева стержня положить поток F(t)=0, то будет происходить
# остывание, пока температура не выровняется по всей длине и не станет равной T0
def F0(t):
    return Fmax * t * exp(-(t / tmax - 1)) / tmax














# Краевые условия
# При х = 0
# def left_boundary_condition(y_prev):
#     y_prev_0 = y_prev[0]
#
#     c_one_half = func_plus_half(c, y_prev_0, tau)
#     c0 = c(y_prev_0)
#     kappa_one_half = func_plus_half(k, y_prev_0, tau)
#     p_one_half =
#     # k_plus = approc_plus_half(k, y_prev_0, tau)
#
#
#     K0 = h / 8 * c_one_half + h / 4 * c0 + t / h * k_plus + \
#          t * h / 8 * p(h / 2) + t * h / 4 * p(0)
#
#     M0 = h / 8 * c_one_half - t / h * k_plus + t * h / 8 * p(h / 2)
#
#     P0 = h / 8 * c_one_half * (y_prev_0 + y_prev[1]) + \
#          h / 4 * c0 * y_prev_0 + F0 * t + t * h / 8 * (3 * f(0) + f(h))
#
#     return K0, M0, P0

EPS = 1e-6

xmin = 0
xmax = l
h = 1e-3

tmin = 0
tmax = tmax
tau = 1e-3


x_array = list(np.arange(xmin, xmax + EPS, h))
t_array = list(np.arange(tmin, tmax + EPS, tau))

N = len(x_array)
M = len(t_array)


# Простая аппроксимация
def func_plus_half(func, x):
    return (func(x) + func(x + h)) / 2

def func_minus_half(func, x, h):
    return (func(x) + func(x - h)) / 2




T = []

# краевое условие в t=0

T.append([T0] * M)
T_0 = [T0] * M


class SimpleIterations:
    # коэффициенты берутся с предыдущего шага
    # то есть первый раз m=1, а считать буду из 0 коэффициентов
    # m = 1,...
    # хз как правильно, пока

    def left_boundary_conditions(self, m):
        def count_K0():
            sum1 = h * func_plus_half(c_T, y0)

    # поступает m, мы будем искать для m+1
    def solve(self, m):
        y_beg_array = T[m]
        return
