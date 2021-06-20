import numpy as np

class pyAirBall:


    # 1. доопределить переменные
    # 2. переделать функции - их логику работы с перемеными
    # 3. где возможно - применить функциоанльное программирование

    def __init__(self, **kwargs):

        pass

    def interface_vector(self, q):
        """
        This function return f vector, that make from q vector
        """
        return np.array([q[0], q[1], q[2] + (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])])


    def interface_u(self):
        """
        This function return an array of gas's speed at interface
        """
        return np.linspace(0, self.vp, self.n + 1)


    def interface_c(self):
        """
        This function return an array of sonic speed at interface
        """
        return (self.c[:self.n + 1] + np.append(self.c[-1], self.c[:self.n])) / 2


    def interface_Max(self, v, u, c):
        """
        This function return an array of Max at interface
        """
        v_left = np.append(v[-1], v[:self.n])
        M_left = (v_left - u) / c
        v_rigth = v[:self.n + 1]
        M_right = (v_rigth - u) / c
        ret = f_plus(M_left) + f_minus(M_right)
        return ret


    def interface_p(self, v, u, c, p):
        """
        This function return an array of pressure at interface
        """
        p_left = np.append(p[-1], p[:n])
        v_left = np.append(v[-1], v[:n])
        M_left = (v_left - u) / c
        v_rigth = v[:n + 1]
        p_rigth = p[:n + 1]
        M_right = (v_rigth - u) / c
        return g_plus(M_left) * p_left + g_minus(M_right) * p_rigth


    def f_plus(self, M):
        """
        This function return an array of Max calculated by f+ scheme
        """
        b = 1 / 8
        return np.where(np.abs(M) >= 1, 0.5 * (M + np.abs(M)), 0.25 * ((M + 1) * (M + 1)) * (1 + 4 * b * (M - 1) * (M - 1)))


    def f_minus(self, M):
        """
        This function return an array of Max calculated by f- scheme
        """
        b = 1 / 8
        return np.where(np.abs(M) >= 1, 0.5 * (M - np.abs(M)),
                        -0.25 * ((M - 1) * (M - 1)) * (1 + 4 * b * (M + 1) * (M + 1)))


    def g_plus(self, M):
        """
        This function return an array of Max calculated by g+ scheme
        """
        a = 3 / 16
        return np.where(np.abs(M) >= 1, (M + np.abs(M)) / (2 * M),
                        ((M + 1) * (M + 1)) * (((2 - M) / 4) + a * M * (M - 1) * (M - 1)))


    def g_minus(self, M):
        """
        This function return an array of Max calculated by g- scheme
        """
        a = 3 / 16
        return np.where(np.abs(M) >= 1, (M - np.abs(M)) / (2 * M),
                        ((M - 1) * (M - 1)) * (((2 + M) / 4) - a * M * (M + 1) * (M + 1)))


    def new_dt(self, dx, v, cq):
        """
        This function return new time step
        """
        dt = ku * np.min(dx / (np.abs(v) + cq))
        return dt


    def p_fun(self, q):
        """
        This function return an array of pressure in each cell
        """
        return (k - 1) * (q[2] - 0.5 * q[1] * q[1] / q[0])


    def c(self, q):
        """
        This function return an array of sonic speed in each cell
        """
        p = p_fun(q)
        return np.sqrt((k / q[0]) * p)
