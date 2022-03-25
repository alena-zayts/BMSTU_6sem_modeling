from math import exp
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

EPS = 1e-4
c = 3e10
k0 = 8e-4
m = 0.786
R = 0.35
Tw = 2000
T0 = 1e4
p_range = [4, 15 + 1]


def T(z, p):
    return (Tw - T0) * (z ** p) + T0


def up(z, p):
    return 3.084e-4 / (exp(4.709e4 / T(z, p)) - 1)


def k(z, p):
    return k0 * ((T(z, p) / 300) ** 2)

def psi(uF_tuple):
    u, F = uF_tuple
    return F - m * c * u / 2


class RungeKutta4Solver:
    def __init__(self, f_func, phi_func):
        self.f = f_func
        self.phi = phi_func

    def next_xyz(self, xn, yn, zn, h):
        k1 = h * self.f(xn, yn, zn)
        q1 = h * self.phi(xn, yn, zn)
        q2 = h * self.phi(xn + h / 2, yn + k1 / 2, zn + q1 / 2)
        k2 = h * self.f(xn + h / 2, yn + k1 / 2, zn + q2 / 2)
        k3 = h * self.f(xn + h / 2, yn + k2 / 2, zn + q2 / 2)
        q3 = h * self.phi(xn + h / 2, yn + k2 / 2, zn + q2 / 2)
        k4 = h * self.f(xn + h, yn + k3, zn + q3)
        q4 = h * self.phi(xn + h, yn + k3, zn + q3)

        x_n_plus_1 = xn + h
        y_n_plus_1 = yn + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        z_n_plus_1 = zn + (q1 + 2 * q2 + 2 * q3 + q4) / 6

        return x_n_plus_1, y_n_plus_1, z_n_plus_1

    def solve(self, x_start, y_start, z_start, x_max, h):
        solution_steps = []
        x, y, z = x_start, y_start, z_start

        while x < x_max + EPS:
            solution_steps.append([x, y, z])
            x, y, z = self.next_xyz(x, y, z, h)

        process_results(solution_steps)

        return tuple((y, z))


class HalfDivisionSolver:
    def __init__(self, f_func):
        self.f = f_func

    def solve(self, x_min, x_max):
        x1 = x_min
        x2 = x_max
        x = (x1 + x2) / 2

        while abs(x1 * x2 / x) > EPS:
            cur_y = self.f(x)
            if abs(cur_y) < EPS:
                return cur_y

            if self.f(x1) * cur_y > 0:
                x1 = x
            else:
                x2 = x

            x = (x1 + x2) / 2

        return x

global P_SHOW

def process_results(solution_steps):
    global P_SHOW
    x, y, z = [step[0] for step in solution_steps], \
              [step[1] for step in solution_steps], \
              [step[2] for step in solution_steps]
    table = pd.DataFrame(index=x)
    table['x (z)'] = x
    table = table.set_index('x (z)')
    table['y (u)'] = y
    table['z (F)'] = z

    round_accuracy = 3
    show_each = 10  # -round_accuracy - 1

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.float_format', lambda x: f'%.{round_accuracy}f' % x)

    print(f'P: {P_SHOW}')
    print(table.iloc[::show_each, :])
    print('\n\n\n')

    fig = plt.figure(figsize=(12, 7))
    plt.xlabel('z')
    # Вывод графиков
    plt.subplot(1, 2, 1)
    plt.plot(x, y, '-', label='u(z)')
    # ups = [up(x0, P_SHOW) for x0 in x]
    ups = list(map(lambda x0: up(x0, P_SHOW), x))
    plt.plot(x, ups, '-', label='up(z)')
    plt.legend()
    plt.subplot(1, 2, 2)
    plt.plot(x, z, '--', label='F(z)')
    plt.savefig(f'p{P_SHOW}.png')
    plt.close(fig)


def main():
    h = 1e-3

    z_min = 1e-7  # zero division
    z_max = 1

    F_min = 0

    ksi_min = 0.01
    ksi_max = 1

    original_stdout = sys.stdout
    os.remove('results.txt')
    with open('results.txt', 'a') as f:
        for p in range(*p_range):
            global P_SHOW
            P_SHOW = p

            sys.stdout = f
            rk_f_func = lambda x, u, v: -3 * R * v * k(x, p) / c
            rk_phi_func = lambda x, u, v: R * c * k(x, p) * (up(x, p) - u) - (v / x)
            u_min_func = lambda ksi: ksi * up(ksi, p)
            hd_f_func = lambda ksi: \
                psi(RungeKutta4Solver(rk_f_func, rk_phi_func).solve(z_min, u_min_func(ksi), F_min, z_max, h))
            solution = HalfDivisionSolver(hd_f_func).solve(ksi_min, ksi_max)

            sys.stdout = original_stdout  #
            print(f'p: {p}, solution: {solution}')


if __name__ == '__main__':
    main()






