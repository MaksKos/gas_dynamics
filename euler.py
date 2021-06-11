
import numpy as np
import time as t

np.seterr(divide='ignore', invalid='ignore')

def interface_vector(q):
    return np.array(
        [q[0],
         q[1],
         q[2] + (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])]
    )


def interface_u(vp):
    return np.linspace(0, vp, n + 1)


def interface_c(c):
    c_right = np.roll(c, 1)
    c_left = c
    c_interf = (c_left + c_right) / 2
    return np.delete(c_interf, n + 1)


def interface_Max(v, u, c):
    v_left = np.roll(v, 1)
    v_left = np.delete(v_left, n + 1)
    M_left = (v_left - u) / c
    v_rigth = np.delete(v, n + 1)
    M_right = (v_rigth - u) / c
    ret = f_plus(M_left) + f_minus(M_right)
    return ret


def interface_p(v, u, c, p):
    p_left = np.roll(p, 1)
    p_left = np.delete(p_left, n + 1)
    v_left = np.roll(v, 1)
    v_left = np.delete(v_left, n + 1)
    M_left = (v_left - u) / c
    v_rigth = np.delete(v, n + 1)
    p_rigth = np.delete(p, n + 1)
    M_right = (v_rigth - u) / c
    return g_plus(M_left) * p_left + g_minus(M_right) * p_rigth


def f_plus(M):
    b = 1 / 8
    return np.where(np.abs(M) >= 1, 0.5 * (M + np.abs(M)), 0.25 * ((M + 1) ** 2) * (1 + 4 * b * (M - 1) ** 2))


def f_minus(M):
    b = 1 / 8
    return np.where(np.abs(M) >= 1, 0.5 * (M - np.abs(M)), -0.25 * ((M - 1) ** 2) * (1 + 4 * b * (M + 1) ** 2))


def g_plus(M):
    a = 3 / 16
    return np.where(np.abs(M) >= 1, (M + np.abs(M)) / (2 * M), ((M + 1) ** 2) * (((2 - M) / 4) + a * M * (M - 1) ** 2))


def g_minus(M):
    a = 3 / 16
    return np.where(np.abs(M) >= 1, (M - np.abs(M)) / (2 * M), ((M - 1) ** 2) * (((2 + M) / 4) - a * M * (M + 1) ** 2))


def new_dt(dx, v, cq):
    dt = ku * np.min(dx / (np.abs(v[:]) + cq[:]))
    return dt


def p_fun(q):
    return (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])


def c(q):
    p = p_fun(q)
    return np.sqrt((k / q[0]) * p)




p0 = 5 *10**6
x0 =0.5
v0 =0
T0 =300
d= 0.03
l = 2
s = np.pi * d ** 2 / 4
m = 0.1
k = 1.4
R = 287
ku = 0.5
n = 300

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

time = 0

t_start = t.time()

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

    F_left=np.roll(F,1,axis = 1)

    F_left = np.delete(F_left,n+1, axis = 1)

    F_right = np.delete(F, n+1, axis = 1)

    f = (c / 2) * (M * (F_right + F_left) - np.abs(M) * (F_right - F_left)) + np.array([np.zeros_like(p), p, p * u])

    q_right = q[: , 0:n + 1]

    q_next = (dx / dx_next) * (q_right - (dt / dx) * (np.roll(f, -1) - f))

    q = q_next

    q = np.append(q, [[0], [0], [0]], axis=1)
    q[0, -1], q[1, -1], q[2, -1] = q[0, 0], -q[1, 0], q[2, 0]
    q[0, n], q[1, n], q[2, n] = q[0, n - 1], -q[1, n - 1] + 2 * vp * q[0, n - 1], q[2, n - 1]

    time += dt
    dx = dx_next

t_end = t.time()
print(t_end - t_start)

print(vp)

print(time)

print(i)








