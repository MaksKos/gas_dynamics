
import numpy as np
import time as t

np.seterr(divide='ignore', invalid='ignore')

def prymay_zadasha(p_function, x0_function):

    def interface_vector(q):
        # nterf_1 = q[0]
        # interf_2 = q[1]
        # interf_3 = q[2] + (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])
        return np.array([q[0], q[1], q[2] + (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])])


    def interface_u(vp):
        return np.linspace(0, vp, n + 1)


    def interface_c(c):
        return (c[:n + 1] + np.append(c[-1], c[:n])) / 2


    def interface_Max(v, u, c):
        v_left = np.append(v[-1], v[:n])
        M_left = (v_left - u) / c
        v_rigth = v[:n + 1]
        M_right = (v_rigth - u) / c
        ret = f_plus(M_left) + f_minus(M_right)
        return ret


    def interface_p(v, u, c, p):
        p_left = np.append(p[-1], p[:n])
        v_left = np.append(v[-1], v[:n])
        M_left = (v_left - u) / c
        v_rigth = v[:n + 1]
        p_rigth = p[:n + 1]
        M_right = (v_rigth - u) / c
        return g_plus(M_left) * p_left + g_minus(M_right) * p_rigth


    def f_plus(M):
        b = 1 / 8
        return np.where(np.abs(M) >= 1, 0.5 * (M + np.abs(M)), 0.25 * ((M + 1) * (M + 1)) * (1 + 4 * b * (M - 1) * (M - 1)))


    def f_minus(M):
        b = 1 / 8
        return np.where(np.abs(M) >= 1, 0.5 * (M - np.abs(M)),
                        -0.25 * ((M - 1) * (M - 1)) * (1 + 4 * b * (M + 1) * (M + 1)))


    def g_plus(M):
        a = 3 / 16
        return np.where(np.abs(M) >= 1, (M + np.abs(M)) / (2 * M),
                        ((M + 1) * (M + 1)) * (((2 - M) / 4) + a * M * (M - 1) * (M - 1)))


    def g_minus(M):
        a = 3 / 16
        return np.where(np.abs(M) >= 1, (M - np.abs(M)) / (2 * M),
                        ((M - 1) * (M - 1)) * (((2 + M) / 4) - a * M * (M + 1) * (M + 1)))


    def new_dt(dx, v, cq):
        dt = ku * np.min(dx / (np.abs(v) + cq))
        return dt


    def p_fun(q):
        return (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])


    def c(q):
        p = p_fun(q)
        return np.sqrt((k / q[0]) * p)


    p0 = p_function
    x0 = x0_function
    v0 = 0
    T0 = 300
    d = 0.0127
    N=40
    l = N*d
    s = np.pi * d ** 2 / 4
    m = 0.011
    k = 1.4
    R = 297
    ku = 0.5
    n = 200

    ro0 = p0 / (R * T0)
    E0 = 1 / (k - 1) * p0 / ro0
    c0 = np.sqrt(k * p0 / ro0)

    q = np.array(3 * [np.zeros(n + 2)], float)
    q[0].fill(ro0)
    q[2].fill(ro0 * E0)
    # Внимани !! n-1 является узлом поршня а n=1 узлом дна
    q[0, -1], q[1, -1], q[2, -1] = q[0, 0], -q[1, 0], q[2, 0]
    q[0, n], q[1, n], q[2, n] = q[0, n - 1], -q[1, n - 1] + 2 * v0 * q[0, n - 1], q[2, n - 1]

    xp = x0
    i = 0
    dx = x0 / n
    vp = 0
    dt = 0
    dx = x0 / n

    #time = 0
   # t_start=t.time()

    while xp < l:
        i += 1

        v = q[1] / q[0]
        pq = p_fun(q)
        cq = np.sqrt((k / q[0]) * pq)

        dt = new_dt(dx, v, cq)

        p_p = p_fun(q)[n - 1]
        vp += (s / m) * dt * p_p
        xp += vp * dt
        dx_next = xp / n

        u = interface_u(vp)
        c = interface_c(cq)
        M = interface_Max(v, u, c)
        p = interface_p(v, u, c, pq)

        F = interface_vector(q)

        F_left = np.append(F[:, n + 1:n + 2], F[:, 0:n + 1], axis=1)

        F_left = F_left[:, :n + 1]

        F_right = F[:, :n + 1]

        f = (c / 2) * (M * (F_right + F_left) - np.abs(M) * (F_right - F_left)) + np.array([np.zeros_like(p), p, p * u])

        q_right = q[:, 0: n + 1]

        q_next = (dx / dx_next) * (q_right - (dt / dx) * (np.roll(f, -1) - f))

        q = q_next

        q = np.append(q, [[0], [0], [0]], axis=1)
        q[0, -1], q[1, -1], q[2, -1] = q[0, 0], -q[1, 0], q[2, 0]
        q[0, n], q[1, n], q[2, n] = q[0, n - 1], -q[1, n - 1] + 2 * vp * q[0, n - 1], q[2, n - 1]

        #time += dt
        dx = dx_next

    #t_end = t.time()
   # print(t_end - t_start)
    return vp

T0 = 300
d = 0.0127
N=40
l = N*d
s = np.pi * d ** 2 / 4
m = 0.011
k = 1.4
R = 297
pmax=10e6

#максимальные и минимальные параметры:
x_min=0
x_max=l/3

p_min=0
p_max=pmax
#количесво значенйи по x0
num_x=10
#количесвто значений по давлению
num_p=10
#время работы прямой задачи
t_prym=1

print("количество точек = ", num_x*num_p)
print("время расчета = ", num_x*num_p*t_prym/60, "минут")


xp=np.linspace(x_min,x_max,num_x)
p0=np.linspace(p_min,p_max,num_p)
m=np.array([])
vp=np.array([])

def m_func (x,p):
    return (s*x*p)/(R*T0)

t_s=t.time()

for i in range(0,num_x):
    print(100*i/num_x, " % ")
    for j in range(0,num_p):
        m=np.append(m,m_func(xp[i],p0[j]))
        vp=np.append(vp,prymay_zadasha(p0[i],xp[i]))


t_end=t.time()

print ('time work:',t_end - t_s)
