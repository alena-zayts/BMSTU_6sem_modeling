from matplotlib import pyplot as plt

c = 3 * 10 ** 10
k_0 = 8 * 1e-3
m = 0.786
r = 0.35
t_w = 2000
t_0 = 10000
p = 4
fig_num = 0
eps = 1e-4

from math import exp

def t(z):
    return (t_w - t_0) * (z ** p) + t_0

def k(z):
    tz = t(z)
    return k_0 * tz * tz / 90000

def u_p(z):
    return 3.084 * 1e-4 / (exp(47990 / t(z)) - 1)


def graph(x, y, title, color):
    global fig_num
    fig_num += 1
    plt.figure(fig_num)
    plt.plot(x, y, color=color)
    # plt.legend()
    plt.title(title)
    plt.grid()

def function_to_file(x, y, name):
    fout = open(name, "w")
    for i in range(0, len(x), 10000):
        fout.write("{:>15.7g}{:>15.7g}\n".format(x[i], y[i]))
    fout.close()

def u_p_show():
    z = 0
    h = 1e-4
    z_s = [0]
    u_s = [u_p(z)]
    u = u_p(z)
    # global fig_num
    # fig_num -= 1
    while z < 1:
        z += h
        z_s.append(z)
        u_s.append(u_p(z))
    graph(z_s, u_s, "Распределение равновесной плотности излучения U_p от z", 'b')
    function_to_file(z_s, u_s, "U-p(x) table")

def kappa(z):
    h = 1e-6
    return c * (k(z + h/2) + k(z - h/2)) / 6 / r / k(z + h/2) / k(z - h/2)

def solve():
    h = 1e-2
    z_half = h / 2

    M0 = kappa(z_half) * z_half + c * r * h * h / 8 * k(z_half) * z_half
    K0 = -kappa(z_half) * z_half + c * r * h * h / 8 * k(z_half) * z_half
    P0 = c * r * h * h / 4 * k(z_half) * u_p(z_half) * z_half
    xi = [0]
    etha = [0]
    xi.append(-K0 / M0)
    etha.append(P0 / M0)
    for i in range(1, int(1 / h) + 1):
        a_n = kappa(z_half) * z_half
        z_half += h/2
        b_n = a_n + c * k(z_half) * z_half * h * h * r
        f_n = c * k(z_half) * u_p(z_half) * h * h * z_half * r
        z_half += h/2
        c_n = kappa(z_half) * z_half
        b_n += c_n 
        xi.append(c_n / (b_n - a_n * xi[i]))
        etha.append((f_n + a_n * etha[i])/(b_n - a_n * xi[i]))

    y = [0 for i in range(len(xi))]
    z_half = 1 - h/2
    MN = -kappa(z_half) * z_half + r * c * h * h / 8 * k(z_half) * z_half
    KN = kappa(z_half) * z_half + m * c * h / 2 + r * c * h * h / 8 * k(z_half) * z_half + r * c * h * h / 4 * k(1)
    PN = c * r * h * h / 4 * (k(z_half) * u_p(z_half) * z_half + k(1) * u_p(1))
    y[-1] = (PN - MN * etha[-1]) / (MN * xi[-1] + KN)
    for i in range(len(y) - 2, -1, -1):
        y[i] = y[i+1] * xi[i+1] + etha[i+1]
    z = [(0 + i * h) for i in range(len(y))]
    function_to_file(z, y, "U(x) table")
    # graph(z, y, "Распределение плотности излучения U от z", 'r')
    return y, z

def get_f_derivatives(y, z):
    h = 1e-6
    dy = []
    dy.append((-3 * y[0] + 4 * y[1] - y[2]) / 2 / h)
    for i in range(1, len(y) - 1):
        dy.append((y[i+1] - y[i-1]) / 2 / h)
    dy.append((3 * y[-1] - 4 * y[-2] + y[-3]) / 2 / h)
    f = [-c / 3 / r / k(z[i]) * dy[i] for i in range(len(z))]
    f[0] = 0

    function_to_file(z, f, "F(x) table")
    graph(z, f, "Распределение плотности излучения F от z", 'r')

    return f

def get_f_integrals(u, z):
    f = [0 for i in range(len(z))]
    h = 1e-6
    s = 0
    z_this = 0
    z_next = h
    for i in range(1, len(f)):
        s += h / 2 * (k(z_this) * (u_p(z_this) - u[i-1]) * z_this + k(z_next) * (u_p(z_next) - u[i]) * z_next)
        f[i] = c * r / z_next * s
        z_this = z_next
        z_next += h

    function_to_file(z, f, "F1(x) table")
    graph(z, f, "Распределение плотности излучения F1 от z", 'r')

    return f


u, z = solve()
f = get_f_derivatives(u, z)
f1 = get_f_integrals(u, z)
# u_p_show()



# x = []
# y = []
# for i in range(5, 16):
#     p = i
#     xi = shoot()
#     u = xi * u_p(0)
#     x.append(p)
#     y.append(u)

# p = 4

# graph(x, y, "Зависимость плотности энергии излучения U(0) от разных значений p")

# x = []
# y = []
# for i in range(5, 505, 50):
#     k_0 = i / 10000
#     xi = shoot()
#     u = xi * u_p(0)
#     x.append(k_0)
#     y.append(u)
# graph(x, y, "Зависимость плотности энергии излучения U(0) от разных значений k_0")
# k_0 = 0.0008


# x = []
# y = []
# for i in range(10, 75, 5):
#     r = i / 100
#     xi = shoot()
#     u = xi * u_p(0)
#     x.append(r)
#     y.append(u)
# graph(x, y, "Зависимость плотности энергии излучения U(0) от разных значений r")
# r = 0.35


# x = []
# y = []
# for i in range(5000, 16000, 1000):
#     t_0 = i
#     xi = shoot()
#     u = xi * u_p(0)
#     x.append(t_0)
#     y.append(u)
# graph(x, y, "Зависимость плотности энергии излучения U(0) от разных значений T_0")
# t_0 = 10000

# x = []
# y = []
# for i in range(500, 6000, 500):
#     t_w = i
#     xi = shoot()
#     u = xi * u_p(0)
#     x.append(t_w)
#     y.append(u)
# graph(x, y, "Зависимость плотности энергии излучения U(0) от разных значений T_w")
# t_w = 2000


plt.show()