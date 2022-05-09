from math import exp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import time

DEBUG = False

# -------------------------------------------------- Constants
a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1


c2 = 0.528e5
m2 = 1

k0 = 1.0
T0 = 300  # температура обдувающего воздуха
R = 0.5  # радиус стержня


# Fmax = 50  # амплитуда импульса (определяет крутизну фронтов)
# tmax = 60  # время (определяет длительность импульса)




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





class Solver:
    # a2 = 2.049
    # b2 = 0.563e-3
    def __init__(self, xmin=0, xmax=10, h=1e-2, tmin=0, tau=1, EPS=None, Fmax=50, tmax=60, a2=2.049, b2=0.563e-3):
        self.xmin = xmin
        self.xmax = xmax
        self.h = h
        self.tmin = tmin
        self.tau = tau
        if EPS:
            self.EPS = EPS
        else:
            self.EPS = self.h / 100
        self.Fmax = Fmax
        self.tmax = tmax
        self.l = xmax - xmin
        self.a2 = a2
        self.b2 = b2

    # коэффициент теплопроводимости материала стержня
    def lambda_T(self, n, m):
        return a1 * (b1 + c1 * pow(self.T[m][n], m1))

    def kappa_T(self, n, m):
        return self.lambda_T(n, m)

    # коэффициент теплоемкости материала стержня
    def c_T(self, n, m):
        return self.a2 + self.b2 * pow(self.T[m][n], m2) - c2 / (self.T[m][n] ** 2)

    # коэффициент оптического поглощения материала стержня
    def k_T(self, n, m):
        return k0 * pow(self.T[m][n] / 300, 2)

    def exp_Tx(self, n, m):
        return exp(-self.k_T(n, m) * self.x_array[n])

    # коэффициент теплоотдачи при обдуве
    def alpha_x(self, n):
        alpha0 = 0.05
        alphaN = 0.01
        d_for_alpha = (alphaN * self.l) / (alphaN - alpha0)
        c_for_alpha = - alpha0 * d_for_alpha
        return c_for_alpha / (self.x_array[n] - d_for_alpha)

    # !?
    # импульс излучения
    # Если в настоящей работе задать поток постоянным, т.е. F0(t) =const, то будет
    # происходить формирование температурного поля от начальной температуры T0
    # до некоторого установившегося (стационарного) распределения T(x,t)
    # Это поле в дальнейшем с течением времени меняться не будет. Это полезный факт для тестирования программы.
    # Если после разогрева стержня положить поток F(t)=0, то будет происходить
    # остывание, пока температура не выровняется по всей длине и не станет равной T0
    def F0_t(self, m):
        t = self.tmin + self.tau * m
        return self.Fmax * t * exp(-(t / self.tmax - 1)) / self.tmax
        # return Fmax * t_array[m] * exp(-(t_array[m] / tmax - 1)) / tmax

    def solve(self):
        self.x_array = list(np.arange(self.xmin, self.xmax + self.EPS, self.h))
        self.N = len(self.x_array)

        self.T = []
        self.t_array = []
        self.T.append([T0] * self.N)
        self.t_array.append(0)

        m = 1
        while True:
            if DEBUG:
                print(f"\n\nStart counting for time {m} ({self.t_array[-1] + self.tau})")
            #print(f"\n\nStart counting for time {m} ({self.t_array[-1] + self.tau})")
            # начальное приближение - сошедшееся с предыдущего шага
            self.T.append(self.T[-1])

            iterations_counter = 0
            while True:
                iterations_counter += 1
                if DEBUG:
                    print(f"Iteration {iterations_counter}")
                if iterations_counter > 1:
                    print('wow, iteration > 1', iterations_counter)
                T_m_prev_iteration = self.T[-1]
                T_m_cur_iteration = self.count_T_array_for_time_m(m)
                self.T[-1] = T_m_cur_iteration
                if self.stop_iterating(T_m_prev_iteration, T_m_cur_iteration):
                    break

            if self.stop_timing():
                break

            m += 1
            self.t_array.append(self.t_array[-1] + self.tau)
        self.t_array.append(self.t_array[-1] + self.tau)
        return self.get_results()


    def get_results(self):
        return self.T, self.t_array

    def check_T(self):
        for Tm in self.T:
            for T_cur in Tm:
                if T_cur > 2000:
                    print('\n\n\n\n\n\n!!!!!!!!!!T>2000\n\n\n\n\n\n\n')


    def count_boundary_conditions(self, m):
        # считаем для того же m, но с предыдущей итерации
        # изначально ошиблась, думая,что с предыдущего времени
        m += 1

        # N должна быть просто последний элемент массива
        self.N -= 1

        K0 = (self.h * func_plus_half(self.c_T, 0, m - 1) / 8 +
              self.h * self.c_T(0, m - 1) / 4 +
              self.tau * func_plus_half(self.kappa_T, 0, m-1) / self.h +
              self.h * self.tau * func_plus_half(self.alpha_x, 0) / (4 * R) +
              self.h * self.tau * self.alpha_x(0) / (2 * R) +
              self.alpha_x(0) * self.tau)

        M0 = (self.h * func_plus_half(self.c_T, 0, m - 1) / 8 -
              self.tau * func_plus_half(self.kappa_T, 0, m-1) / self.h +
              self.h * self.tau * func_plus_half(self.alpha_x, 0) / (4 * R))

        P0 = (self.h * func_plus_half(self.c_T, 0, m - 1) * (self.T[m - 1][0] + self.T[m - 1][1]) / 8 +
              self.h * self.c_T(0, m - 1) * self.T[m - 1][0] / 4 +
              self.alpha_x(0) * self.tau * T0 +
              self.h * self.tau * (
                    self.F0_t(m - 1) * (self.k_T(0, m - 1) * self.exp_Tx(0, m - 1) + func_plus_half(self.k_T, 0, m - 1) * func_plus_half(self.exp_Tx, 0, m - 1)) +
                    2 * T0 * (func_plus_half(self.alpha_x, 0) + self.alpha_x(0)) / R
              ) / 4)

        KN = (self.h * func_minus_half(self.c_T, self.N, m - 1) / 8 +
              self.h * self.c_T(self.N, m - 1) / 4 +
              self.tau * func_minus_half(self.kappa_T, self.N, m-1) / self.h +
              self.h * self.tau * func_minus_half(self.alpha_x, self.N) / (4 * R) +
              self.h * self.tau * self.alpha_x(self.N) / (2 * R) +
              self.alpha_x(self.N) * self.tau)


        MN = (self.h * func_minus_half(self.c_T, self.N, m - 1) / 8 -
              self.tau * func_minus_half(self.kappa_T, self.N, m-1) / self.h +
              self.h * self.tau * func_minus_half(self.alpha_x, self.N) / (4 * R))

        PN = (self.h * func_minus_half(self.c_T, self.N, m - 1) * (self.T[m - 1][self.N] + self.T[m - 1][self.N - 1]) / 8 +
              self.h * self.c_T(self.N, m - 1) * self.T[m - 1][self.N] / 4 +
              self.alpha_x(self.N) * self.tau * T0 +
              self.h * self.tau * (
                    self.F0_t(m - 1) * (self.k_T(self.N, m - 1) * self.exp_Tx(self.N, m - 1) + func_minus_half(self.k_T, self.N, m - 1) * func_minus_half(self.exp_Tx, self.N, m - 1)) +
                    2 * T0 * (func_minus_half(self.alpha_x, self.N) + self.alpha_x(self.N)) / R
              ) / 4)

        self.N += 1
        return {"K0": K0, "M0": M0, "P0": P0,
                "KN": KN, "MN": MN, "PN": PN,}

    def count_coeffs(self, n, m):
        # считаем для того же m, но с предыдущей итерации
        # изначально ошиблась, думая,что с предыдущего времени
        m += 1

        A = self.tau * func_minus_half(self.kappa_T, n, m - 1) / self.h
        D = self.tau * func_plus_half(self.kappa_T, n, m - 1) / self.h
        B = self.h * self.c_T(n, m - 1) + A + D + 2 * self.h * self.tau * self.alpha_x(n) / R
        F = self.h * self.c_T(n, m - 1) * self.T[m - 1][n] + self.h * self.tau * self.k_T(n, m - 1) * self.F0_t(m - 1) * self.exp_Tx(n, m - 1) + 2 * self.h * self.tau * T0 * self.alpha_x(n) / R
        return {"A": A, "B": B, "D": D, "F": F}

    def count_T_array_for_time_m(self, m):
        boundary_conditions = self.count_boundary_conditions(m)

        # Прямой ход
        ksi_array = [-boundary_conditions['M0'] / boundary_conditions["K0"]]
        eta_array = [boundary_conditions['P0'] / boundary_conditions["K0"]]

        # x = self.h
        # n = 1
        for n in range(1, self.N - 1):
        #while (x + self.h < l):
            cur_coeffs = self.count_coeffs(n, m)
            ksi_array.append(cur_coeffs["D"] / (cur_coeffs["B"] - cur_coeffs["A"] * ksi_array[n - 1]))
            eta_array.append((cur_coeffs["F"] + cur_coeffs["A"] * eta_array[n - 1]) /
                             (cur_coeffs["B"] - cur_coeffs["A"] * ksi_array[n - 1]))
            # n += 1
            # x += self.h

        # Обратный ход
        T_cur_array = [0] * self.N
        T_cur_array[-1] = ((boundary_conditions['PN'] - boundary_conditions['MN'] * eta_array[-1]) /
                 (boundary_conditions['KN'] + boundary_conditions['MN'] * ksi_array[-1]))

        for i in range(self.N - 2, -1, -1):
            T_cur_array[i] = ksi_array[i] * T_cur_array[i + 1] + eta_array[i]

        return T_cur_array

    def stop_timing(self):
        T_array_prev, T_array_new = self.T[-2], self.T[-1]
        for T_prev, T_cur in zip(T_array_prev, T_array_new):
            if abs((T_prev - T_cur) / T_cur) > 10e-4:
                return False
        return True

    def stop_iterating(self, T_array_prev, T_array_new):
        if len(T_array_prev) != self.N or len(T_array_new) != self.N:
            raise ValueError('NNNN')
        differences = [abs((T_array_new[n] - T_array_prev[n]) / T_array_new[n]) for n in range(self.N)]
        max_difference = max(differences)
        return max_difference < 1




    def print_results(self):
        print(self.t_array)
        print(self.x_array)
        for row in self.T:
            print(row)

        columns = [f'{t}' for t in self.t_array]
        table = pd.DataFrame(columns=columns, index=self.x_array)
        for i in range(len(self.t_array)):
            table[f'{self.t_array[i]}'] = self.T[i]

        # table = table.set_index('z')

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_colwidth', None)
        pd.options.display.float_format = '{:,.9e}'.format
        show_each = 1
        print(table.iloc[::show_each, :])
        print('\n\n\n')
        table.to_excel('result.xls')

    def show_results(self):
        M = len(self.T)# обрезаем неинтересную часть графика

        # Трехмерный
        x, y = np.mgrid[0:self.t_array[-1]+self.EPS:self.tau, 0:self.xmax+self.EPS:self.h]
        z = np.array([np.array(T_m) for T_m in self.T])

        fig_3d = plt.figure()
        xyz = fig_3d.add_subplot(111, projection='3d')
        xyz.plot_surface(x, y, z, cmap='Oranges')
        plt.xlabel("t, c")
        plt.ylabel("x, см")
        fig_3d.show()
        plt.savefig(f'3d.png')

        # Проекции
        fig, (first_graph, second_graph) = plt.subplots(
            nrows=1, ncols=2,
            figsize=(8, 4))

        # Первая cm - K
        step_by_time = int(len(self.T) / 10)
        for T_m in self.T[len(self.T) - 1:0:-step_by_time]:
            first_graph.plot(self.x_array, T_m)
        first_graph.set_xlabel("x, cm")
        first_graph.set_ylabel("T, K")
        first_graph.grid()

        # Вторая sec - K
        step_by_cm = int(len(self.x_array) / 10)
        for i in range(0, len(self.x_array), step_by_cm):
            line = [T_m[i] for T_m in self.T]
            second_graph.plot(self.t_array, line)

        second_graph.set_xlabel("t, sec")
        second_graph.set_ylabel("T, K")
        second_graph.grid()
        fig.show()
        plt.savefig(f'2d.png')
        print()

    def check_steps_ok(self):
        eps = 1e-2
        M = len(self.T)
        F0_value = self.alpha_x(0) * (T0 - self.T[-1][0])
        FN_value = self.alpha_x(self.N - 1) * (self.T[-1][-1] - T0)

        # средние
        upper_integral = 0
        for i in range(1, self.N):
            sum1 = self.alpha_x(i) * (self.T[-1][i] - T0)
            sum2 = self.alpha_x(i - 1) * (self.T[-1][i - 1] - T0)
            upper_integral += (sum1 + sum2)
            # upper_integral += func_minus_half(self.alpha_x, i) * (((self.T[-1][i] + self.T[-1][i-1]) / 2) - T0)
        upper_integral *= (self.h / R)
        # upper_integral_2 = scipy.integrate.trapezoid

        F0t = self.F0_t(M - 1)

        lower_integral = 0
        for i in range(1, self.N):
            sum1 = self.k_T(i, M - 1) * self.exp_Tx(i, M - 1)
            sum2 = self.k_T(i - 1, M - 1) * self.exp_Tx(i - 1, M - 1)
            lower_integral += (sum1 + sum2)
            # lower_integral += func_minus_half(self.k_T, i, M - 1) * func_minus_half(self.exp_Tx, i, M-1) * self.h

        chisl = -F0_value + FN_value + upper_integral
        znam = F0t * lower_integral

        counted = abs((chisl / znam) - 1)
        if znam < 1e-6:
            print()
        #  print('counted', counted)
        return counted <= eps, counted



def study_steps():
    print('research')
    h_array = [1, 0.5, 0.1, 0.01, 0.001]
    tau_array = [1, 0.5, 0.1]
    results = pd.DataFrame(columns=['h', 'tau', 'balanced', 'M', 'counted'])
    for h in h_array:
        for tau in tau_array:
            solver = Solver(h=h, tau=tau)
            solver.solve()
            T, t_array = solver.get_results()
            balanced, counted = solver.check_steps_ok()
            print(f'h={h}, tau={tau}, balanced={balanced}, M={len(t_array)}, counted={counted}')
            results = results.append({'h':h, 'tau': tau, 'balanced': balanced, 'M': len(t_array), 'counted': counted}, ignore_index=True)
    print(results)
    results.to_excel('study_steps.xls')


# График зависимости температуры T(0,t) от a2, b2
# С ростом теплоемкости (этих коэффициентов) темп нарастания температуры снижается
# коэффициент теплоемкости материала стержня
# def c_T(self, n, m):
#       return self.a2 + self.b2 * pow(self.T[m][n], m2) - c2 / (self.T[m][n] ** 2)
def study_ab():
    print('research ab')
    a2_array = [2.049, 3, 5]
    b2_array = [0.563e-3, 0.6e-2, 0.6e-1]
    q=0
    fig = plt.figure(figsize=(12, 7))


    for i, a2 in enumerate(a2_array):
        for j, b2 in enumerate(b2_array):
            plt.subplot(len(a2_array), len(b2_array), i * len(b2_array) + j + 1)
            q+=1
            solver = Solver(a2=a2, b2=b2)
            solver.solve()
            T, t_array = solver.get_results()

            T_0t = [T_m[0] for T_m in T]
            plt.plot(t_array, T_0t)

            print(f'a2={a2}, b2={b2}, T_0t={T_0t}')
            title = f"{q} a2={a2} b2={b2}"
            plt.title(title)

    plt.savefig(f"a2b2.png")
    plt.show()






def main():
    # solver = Solver(h=0.0001)
    # solver.solve()
    # solver.show_results()
    # #solver.print_results()
    # print(solver.check_steps_ok())

    #study_steps()
    study_ab()



if __name__ == '__main__':
    main()



