from math import exp
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
import numpy as np

EPS = 1e-4
c = 3e10
k0 = 8e-3  # k0 = 0.01 # k0 = 0.008
m = 0.786
R = 0.35
Tw = 2000
T0 = 1e4
# p_range = [4, 15 + 1]
p = 4


def T(z):
    return (Tw - T0) * (z ** p) + T0


def up(z):
    if z > 1:
        print(z, 3.084e-4 / (exp(4.799e4 / T(z)) - 1))
    return 3.084e-4 / (exp(4.799e4 / T(z)) - 1)


def k(z):
    return k0 * ((T(z) / 300) ** 2)


class Solver:
    def __init__(self, k_func, p_func, f_func):
        self.k_func = k_func
        self.p_func = p_func
        self.f_func = f_func

    # alpha = mc / 2, betta = 0
    def solve(self, z_min, z_max, h, F0=0, alpha=(m * c / 2), beta=0):
        self.h = h
        self.z_array = list(np.arange(z_min, z_max + EPS, h))
        self.zF_array = [z_min] + list(np.arange(z_min + h / 2, z_max, h)) + [z_max]
        self.N = len(self.z_array)# -1
        self.z_array.append(z_max + h)
        self.F0 = F0
        self.alpha = alpha
        self.beta = beta

        self.count_coeffs()
        self.forward()
        self.bacward()
        self.count_F_array()
        return self.z_array[:-1], self.zF_array, self.y_array, self.F_array

    def forward(self):
        self.ksi_array = []
        self.eta_array = []

        #self.ksi_array.append(-self.M0 / self.K0) - my
        ## !!!!!!!
        # self.ksi_array.append(-self.M0 / self.K0)
        # self.eta_array.append(self.P0 / self.K0)
        # M и K перепутала
        # round_to = 6
        self.ksi_array.append(-self.K0 / self.M0)
        self.eta_array.append(self.P0 / self.M0)

        for i in range(1, self.N):
            # self.ksi_array[-1] = round(self.ksi_array[-1], round_to)
            # self.eta_array[-1] = round(self.eta_array[-1], round_to)
            self.ksi_array.append(self.ksi_n_plus_1(i))
            self.eta_array.append(self.eta_n_plus_1(i))
        # self.ksi_array[-1] = round(self.ksi_array[-1], round_to)
        # self.eta_array[-1] = round(self.eta_array[-1], round_to)

        # self.ksi_array[-1] = 0
        # self.eta_array[-1] = 1e-6

        print(f'ksi {len(self.ksi_array)}: {self.ksi_array}')
        print(f'eta {len(self.eta_array)}: {self.eta_array}')
        print(self.N)


    def bacward(self):
        self.y_array = [0] * self.N
        # M и K перепутала
        # self.y_array[-1] = (self.PN - self.KN * self.eta_array[-1]) / (self.KN * self.ksi_array[-1] + self.MN) - не те mn
        # self.y_array[-1] = (self.PN - self.MN * self.eta_array[-1]) / (self.MN * self.ksi_array[-1] + self.KN) - наверное правильно
        # self.y_array[-1] = (self.PN - self.MN * 1e-6) / (self.MN * 0 + self.KN)


        # self.y_array[-1] = (self.PN - self.MN * self.eta_array[-1]) / (self.MN * self.ksi_array[-1] + self.KN)
        # self.y_array[-1] = (self.PN - self.MN * 1e-6) / (self.MN * 0 + self.KN)
        self.y_array[-1] = (self.PN - self.KN * self.eta_array[-1]) / (self.KN * self.ksi_array[-1] + self.MN)
        # x1 = (self.PN - self.KN * self.eta_array[-1]) / (self.KN * self.ksi_array[-1] + self.MN)
        # x2 = (self.PN - self.MN * self.eta_array[-1]) / (self.MN * self.ksi_array[-1] + self.KN)

        for i in range(self.N - 1 - 1, -1, -1):
            # self.y_array[i] = self.ksi_array[i + 1] * self.y_array[i + 1] + self.eta_array[i + 1]
            self.y_array[i] = self.ksi_array[i] * self.y_array[i + 1] + self.eta_array[i]
            # print(f'i={i}, y[i]= {self.ksi_array[i]} * {self.y_array[i + 1]} + {self.eta_array[i]}')

    def count_coeffs(self):
        def count_M0():
            sum1 = self.z_n_plus_half(0) * self.kappa_n_plus_half(0) / (R * R * self.h)
            sum2 = self.h * self.p_n_plus_half(0) * self.V_n_one_halph() / 8
            sum3 = self.h * self.p_func(self.z_array[0]) * self.V_n_one_halph() / 4
            return -(sum1 + sum2 + sum3)

        def count_K0():
            sum1 = self.z_n_plus_half(0) * self.kappa_n_plus_half(0) / (R * R * self.h)
            sum2 = self.h * self.p_n_plus_half(0) * self.V_n_one_halph() / 8
            return sum1 - sum2

        def count_P0():
            sum1 = self.z_array[0] * self.F0 / R
            sum2 = self.h * self.V_n_one_halph() * (self.f_func(self.z_array[0]) + self.f_n_plus_half(0)) / 4
            return -(sum1 + sum2)

        def count_KN():
            sum1 = self.z_n_minus_half(self.N) * self.kappa_n_minus_half(self.N) / (R * R * self.h)
            sum2 = self.h * self.p_n_minus_half(self.N) * self.V_N_minus_half() / 8
            return (sum1 - sum2)

        def count_MN():
            sum1 = self.z_n_minus_half(self.N) * self.kappa_n_minus_half(self.N) / (R * R * self.h)
            sum2 = self.z_array[self.N] * self.alpha / R
            sum3 = self.h * self.p_n_minus_half(self.N) * self.V_N_minus_half() / 8
            sum4 = self.h * self.p_func(self.z_array[self.N]) * self.V_N_minus_half() / 4
            return -(sum1 + sum2 + sum3 + sum4)

        def count_PN():
            sum1 = self.z_array[self.N] * self.alpha * self.beta / R
            # sum1 = 0  # betta = 0
            sum2 = (self.h * self.V_N_minus_half() *
                    (self.f_func(self.z_array[self.N]) + self.f_n_minus_half(self.N)) / 4)
            return -(sum1 + sum2)

        self.M0 = count_M0()
        self.K0 = count_K0()
        self.P0 = count_P0()

        # self.N -= 1
        self.MN = count_MN()
        self.KN = count_KN()
        self.PN = count_PN()
        # self.N += 1

    def count_F_array(self):
        self.F_array = [self.F0]
        for n in range(self.N - 1):
            F_plus_half = self.kappa_n_plus_half(n) * (self.y_array[n] - self.y_array[n + 1]) / (R * self.h) # ?R
            self.F_array.append(F_plus_half)

        self.F_array.append(self.alpha * (self.y_array[-1] - self.beta))

    def func_n_minus_half(self, func, n):
        return float(func(self.z_array[n - 1]) + func(self.z_array[n])) / 2

    def func_n_plus_half(self, func, n):
        return self.func_n_minus_half(func, n + 1)

    def z_n_minus_half(self, n):
        return self.func_n_minus_half(lambda x: x, n)

    def z_n_plus_half(self, n):
        return self.func_n_plus_half(lambda x: x, n)

    def kappa_n_minus_half(self, n):
        return self.func_n_minus_half(self.k_func, n)

    def kappa_n_plus_half(self, n):
        return self.func_n_plus_half(self.k_func, n)

    def p_n_minus_half(self, n):
        return self.func_n_minus_half(self.p_func, n)

    def p_n_plus_half(self, n):
        return self.func_n_plus_half(self.p_func, n)

    def f_n_minus_half(self, n):
        return self.func_n_minus_half(self.f_func, n)

    def f_n_plus_half(self, n):
        return self.func_n_plus_half(self.f_func, n)



    def V_n_one_halph(self):
        return (self.z_n_plus_half(0) ** 2 - self.z_array[0] ** 2) / 2

    def V_N_minus_half(self):
        return (self.z_array[self.N] ** 2 - self.z_n_minus_half(self.N) ** 2) / 2


    def Vn(self, n):
        return (self.z_n_plus_half(n) ** 2 - self.z_n_minus_half(n) ** 2) / 2

    def An(self, n):
        return (self.z_n_minus_half(n) * self.kappa_n_minus_half(n)) / (R * R * self.h)

    def Cn(self, n):
        return (self.z_n_plus_half(n) * self.kappa_n_plus_half(n)) / (R * R * self.h)

    def Bn(self, n):
        return self.An(n) + self.Cn(n) + self.p_func(self.z_array[n]) * self.Vn(n)

    def Dn(self, n):
        return self.f_func(self.z_array[n]) * self.Vn(n)




    def ksi_n_plus_1(self, n):
        return self.Cn(n) / (self.Bn(n) - self.An(n) * self.ksi_array[n - 1])

    def eta_n_plus_1(self, n):
        return (self.Dn(n) + self.An(n) * self.eta_array[n - 1]) / (self.Bn(n) - self.An(n) * self.ksi_array[n - 1])



def process_results(z, zF, u, F):
    print('z    u')
    for z0, u0 in zip(z, u):
        print(z0, u0)

    print('z    F')
    for z0, F0 in zip(zF, F):
        print(z0, F0)

    # F.append('-')
    # table = pd.DataFrame(index=z)
    # table['z'] = z
    # table = table.set_index('z')
    # table['u'] = u
    # table['F'] = F

    # show_each = 1
    #
    # pd.set_option('display.max_rows', None)
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_colwidth', None)
    #
    #
    # print(table.iloc[::show_each, :])
    # print('\n\n\n')

    fig = plt.figure(figsize=(12, 7))

    plt.subplot(1, 2, 1)
    plt.plot(z, u, '-', label='u(z)')

    ups = list(map(lambda z0: up(z0), z))
    plt.plot(z, ups, '-', label='up(z)')
    plt.legend()
    plt.xlabel('z')
    plt.grid()

    plt.subplot(1, 2, 2)
    # w, e = z_reformatted(z), F[:-1]
    # print(w)
    # print(e)
    plt.plot(zF, F, '--', label='F(z)')
    plt.legend()
    plt.xlabel('z')
    plt.grid()
    plt.savefig(f'graph.png')
    plt.show()
    plt.close(fig)


# !!!!
def k_func(z):
    return c / (3 * k(z) * R)

def p_func(z):
    return c * k(z)

def f_func(z):
    return c * k(z) * up(z)


def z_reformatted(z):
    new_z = []
    for i in range(len(z) - 1):
        new_z.append((z[i] + z[i + 1]) / 2)
    return new_z


def main():
    h = 1e-4

    z_min = 0
    z_max = 1

    solver = Solver(k_func, p_func, f_func)
    z, zF, u, F = solver.solve(z_min, z_max, h)
    process_results(z, zF, u, F)




# def main():
#     h = 1e-1
#
#     z_min = 0
#     z_max = 1
#
#     F_min = 0
#
#     ksi_min = 0.01
#     ksi_max = 1
#
#     original_stdout = sys.stdout
#     os.remove('results.txt')
#     with open('results.txt', 'a') as f:
#         for p in range(*p_range):
#             rk_f_func = lambda x, u, v: -3 * R * v * k(x, p) / c
#             rk_phi_func = lambda x, u, v: (
#                 R * c * k(x, p) * (up(x, p) - u) - (v / x) if x != 0  # zero division
#                 else
#                 R * c * k(x, p) * (up(x, p) - u) / 2)
#             u_min_func = lambda ksi: ksi * up(0, p)
#
#             hd_f_func = lambda ksi: \
#                 psi(RungeKutta4Solver(rk_f_func, rk_phi_func).solve(z_min, u_min_func(ksi), F_min, z_max, h))
#
#             ksi = HalfDivisionSolver(hd_f_func).solve(ksi_min, ksi_max)
#
#             global P_SHOW
#             P_SHOW = p
#             sys.stdout = f
#             RungeKutta4Solver(rk_f_func, rk_phi_func).solve(z_min, u_min_func(ksi), F_min, z_max, h, show_result=True)
#             sys.stdout = original_stdout  #
#             print(f'p: {p}, solution_ksi: {ksi}')


if __name__ == '__main__':
    main()