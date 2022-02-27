from prettytable import PrettyTable
import pandas as pd
import matplotlib.pyplot as plt
EPS = 1e-6


class UDESolver:
    def __init__(self, x_start, y_start, x_max, step, f, f_derivatives):
        if (x_start < x_max) != (step > 0):
            raise ValueError('Ошибка в шаге')

        self.x_start = x_start
        self.y_start = y_start
        self.x_max = x_max
        self.step = step
        self.f = f
        self.f_derivatives = f_derivatives
        self.cmp_func = lambda x1, x2: x1 < x2 + EPS

    def reverse_move(self):
        self.step *= -1
        self.x_max *= -1
        self.cmp_func = lambda x1, x2: x1 > x2 - EPS

    def x_range(self):
        result = []
        x = self.x_start
        while self.cmp_func(x, self.x_max):
            result.append(x)
            x += self.step
        return result

    def solve_euler(self):
        result = []
        x, y = self.x_start, self.y_start

        while self.cmp_func(x, self.x_max):
            result.append(y)

            y = y + self.step * self.f(x, y)
            x += self.step

        return result

    def solve_runge_kutta(self):
        a = 0.5
        result = []
        x, y = self.x_start, self.y_start

        while self.cmp_func(x, self.x_max):
            result.append(y)

            k1 = self.f(x, y)
            k2 = self.f(x + self.step / (2 * a), y + self.step * k1 / (2 * a))
            y += self.step * ((1 - a) * k1 + a * k2)
            x += self.step

        return result

    def solve_picar(self, approx):
        func = self.f_derivatives[approx - 1]

        result = []
        x, y = self.x_start, self.y_start

        while self.cmp_func(x, self.x_max):
            result.append(y)
            x += self.step
            y = func(x)

        return result


def function(x, u):
    return x * x + u * u


def fd1(x):
    return pow(x, 3) / 3


def fd2(x):
    return fd1(x) + pow(x, 7) / 63


def fd3(x):
    return fd2(x) + 2 * pow(x, 11) / 2079 + pow(x, 15) / 59535


def fd4(x):
    return (fd2(x) + 2 * pow(x, 11) / 2079 + 13 * pow(x, 15) / 218295 +
            82 * pow(x, 19) / 37328445 + 662 * pow(x, 23) / 10438212015 +
            4 * pow(x, 27) / 3341878155 + pow(x, 31) / 109876903905)

def draw_plots(table):
    plt.figure(figsize=(30, 10))
    x = table.index
    for column_name in table.columns:
        y = table[column_name]
        plt.plot(x, y, label=column_name)

    plt.legend()
    plt.show()

def main():
    x_start = 0
    y_start = 0
    x_max = 1
    step_accuracy = 2
    step = float(f'1e-{step_accuracy}')
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.float_format', lambda x: f'%.{step_accuracy}f' % x)

    solver = UDESolver(x_start, y_start, x_max, step, function, [fd1, fd2, fd3, fd4])

    table = pd.DataFrame(index=solver.x_range())
    table['x'] = solver.x_range()
    table = table.set_index('x')
    table['Euler'] = solver.solve_euler()
    table['Runge-Kutta'] = solver.solve_runge_kutta()
    for i in range(4):
        table[f"Picard, {i + 1}"] = solver.solve_picar(i + 1)


    # print(table)
    # draw_plots(table)


    solver.reverse_move()

    table2 = pd.DataFrame(index=solver.x_range())
    table2['x'] = solver.x_range()
    table2 = table2.set_index('x')
    table2['Euler'] = solver.solve_euler()
    table2['Runge-Kutta'] = solver.solve_runge_kutta()
    for i in range(4):
        table2[f"Picard, {i + 1}"] = solver.solve_picar(i + 1)

    # print(table2)
    # draw_plots(table2)

    full_table = pd.concat([table, table2], sort=True, axis=0)
    full_table = full_table.sort_index(ascending=True)

    print(full_table)
    # draw_plots(full_table)



if __name__ == "__main__":
    main()


# Подбираем шаг:

# Euler:
# При 1e-1 y(1) = 0.2925421046
# При 1e-2 y(1) = 0.3331073593
# При 1e-3 y(1) = 0.3484859823
# При 1e-4 y(1) = 0.3501691515
# При 1e-5 y(1) = 0.3502255745
# Шаг ничего не меняет (между 1e-3 и 1e-4)
# Значит мы подобрали нужный нам шаг.

# Runge:
# При 1e-1 y(1) = 0.3485453439
# При 1e-2 y(1) = 0.3391265967
# При 1e-3 y(1) = 0.3491103993
# При 1e-4 y(1) = 0.3502318426
# При 1e-5 y(1) = 0.3502318443
# Аналогично.
