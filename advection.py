import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm

def ux_exact(x, t):
	return np.exp(2j * np.pi *  (x - c * t))

def dU_dt(U, t):
	ux = U[0:]
	dux_dt = np.zeros(nx, dtype=complex)
	dux_dt[1:-1] = -c * (ux[2:] - ux[0:-2]) / (2 * dx)
	dux_dt[0]    = -c * (ux[1]  - ux[-1]  ) / (2 * dx)
	dux_dt[-1]   = -c * (ux[0]  - ux[-2]  ) / (2 * dx)
	dU_dt = dux_dt
	return dux_dt

def t_integral(U0, t):
	r = complex_ode(dU_dt)
	r.set_initial_value(U0, t_min)


nx = 128
lx = 0.5
x_min = -lx
x_max = lx
dx = (x_max - x_min) / nx
x = np.linspace(x_min + dx / 2, x_max - dx / 2, nx)

nt = 128
t_min = 0
t_max =  1
dt = (t_max - t_min) / (nt - 1)
t = np.linspace(t_min, t_max, nt)

X, T = np.meshgrid(x, t)

c = 1

ux0 = ux_exact(x, t_min)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(X, T, np.real(ux_exact(X,T)), cmap=cm.cool)
# ax.set_xlabel('x')
# ax.set_ylabel('z')


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(X, T, np.imag(ux_exact(X,T)), cmap=cm.cool)
# ax.set_xlabel('x')
# ax.set_ylabel('t')