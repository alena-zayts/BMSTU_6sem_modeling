from prettytable import PrettyTable

EPS = 1e-6


class UDESolver:
    def __init__(self, x_start, y_start, x_max, step, f, f_derivatives):
        self.x_start = x_start
        self.y_start = y_start
        self.x_max = x_max
        self.step = step
        self.f = f
        self.f_derivatives = f_derivatives

    def x_range(self):
        result = []
        x = self.x_start
        while x < self.x_max + EPS:
            result.append(round(x, 2))
            x += self.step
        return result

    def solve_euler(self):
        result = []
        x, y = self.x_start, self.y_start

        while x < self.x_max + EPS:
            result.append(y)

            y = y + self.step * self.f(x, y)
            x += self.step

        return result

    def solve_runge_kutta(self):
        result = []
        x, y = self.x_start, self.y_start

        while x < self.x_max + EPS:
            result.append(y)

            k1 = self.f(x, y)
            k2 = self.f(x + self.step, y + self.step * k1)
            y = y + self.step * (k1 + k2) / 2
            x += self.step

        return result

    def solve_picar(self, approx):
        func = self.f_derivatives[approx - 1]

        result = []
        x, y = self.x_start, self.y_start

        while x < self.x_max + EPS:
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


def main():
    x_start = 0
    y_start = 0
    x_max = 1
    step = 1e-4

    solver = UDESolver(x_start, y_start, x_max, step, function, [fd1, fd2, fd3, fd4])

    tb = PrettyTable()
    tb.add_column("X", solver.x_range())
    tb.add_column("Euler", solver.solve_euler())
    tb.add_column("Runge-Kutta", solver.solve_runge_kutta())
    for i in range(4):
        tb.add_column(f"Picard, {i}", solver.solve_picar(i))

    print(tb)


if __name__ == "__main__":
    main()


# def Euler(x_max, h):
#     result = []
#     x, y = x_start, y_start
#
#     while x < x_max + EPS:
#         result.append(y)
#
#         y = y + h * f(x, y)
#         x += h
#
#     return result
#
#
# def RungeKutta(x_max, h):
#     result = []
#     x, y = x_start, y_start
#
#     while x < x_max:
#         result.append(y)
#
#         k1 = f(x, y)
#         k2 = f(x + h, y + h * k1)
#         y = y + h * (k1 + k2) / 2
#         x += h
#
#     return result
#
# def Picar(x_max, h, depth):
#     if depth == 1:
#         func = fd1
#     elif depth == 2:
#         func = fd2
#     elif depth == 3:
#         func = fd2
#     else:
#         func = fd2
#
#     result = []
#     x, y = x_start, y_start
#
#     while x < x_max:
#         result.append(y)
#         x += h
#         y = func(x)
#
#     return result
#
# def main1():
#     column_names = ["X", "Picard 1", "Picard 2", "Picard 3", "Picard 4", "Runge"]
#
#     tb = PrettyTable()
#     tb.add_column("X", x_range(MAX_X, STEP))
#     tb.add_column("Picard 1", Picar(MAX_X, STEP, 1))
#     tb.add_column("Picard 2", Picar(MAX_X, STEP, 2))
#     tb.add_column("Picard 3", Picar(MAX_X, STEP, 3))
#     tb.add_column("Picard 4", Picar(MAX_X, STEP, 4))
#     tb.add_column("Euler", Euler(MAX_X, STEP))
#     tb.add_column("Runge", RungeKutta(MAX_X, STEP))
#
#     print(tb)


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
