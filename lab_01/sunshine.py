from prettytable import PrettyTable
from decimal import Decimal

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

MAX_X = 1
STEP = 1e-4


def f(x, y):
    return x * x + y * y


# return pow(x, 2) + pow(y, 2)

def fp1(x):
    return pow(x, 3) / 3


def fp2(x):
    return pow(x, 7) / 63 + \
           fp1(x)


def fp3(x):
    return pow(x, 15) / 59535 + \
           2 * pow(x, 11) / 2079 + \
           fp2(x)


def fp4(x):
    return pow(x, 31) / 109876903905 + \
           4 * pow(x, 27) / 3341878155 + \
           662 * pow(x, 23) / 10438212015 + \
           82 * pow(x, 19) / 37328445 + \
           fp3(x)


def Picar(x_max, h, func):
    result = list()
    x, y = 0, 0

    while x < x_max:
        result.append(y)
        x += h
        y = func(x)

    return result


def Euler(x_max, h):
    result = list()
    x, y = 0, 0  # Начальное условие.

    while x < x_max:
        result.append(y)
        # print(y)
        y = y + h * f(x, y)
        x += h

    return result


def Runge(x_max, h):
    result = list()
    coeff = h / 2
    x, y = 0, 0

    while x < x_max:
        result.append(y)
        # print(x, y)
        y = y + h * f(x + coeff, y + coeff * f(x, y))
        x += h

    return result


def x_range(x_max, h):
    result = list()
    x = 0
    while x < x_max:
        result.append(round(x, 2))
        x += h
    return result


def main():
    column_names = ["X", "Picard 1", "Picard 2", "Picard 3", "Picard 4", "Runge"]

    tb = PrettyTable()
    tb.add_column("X", x_range(MAX_X, STEP))
    tb.add_column("Picard 1", Picar(MAX_X, STEP, fp1))
    tb.add_column("Picard 2", Picar(MAX_X, STEP, fp2))
    tb.add_column("Picard 3", Picar(MAX_X, STEP, fp3))
    tb.add_column("Picard 4", Picar(MAX_X, STEP, fp4))
    tb.add_column("Euler", Euler(MAX_X, STEP))
    tb.add_column("Runge", Runge(MAX_X, STEP))

    print(tb)


if __name__ == "__main__":
    main()
