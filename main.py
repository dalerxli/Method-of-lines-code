import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import control as ct

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(ct.X, ct.Z, np.real(ct.ux_ana(ct.X,ct.Z)))
ax.set_xlabel('x')
ax.set_ylabel('z')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(ct.X, ct.Z, np.imag(ct.ux_ana(ct.X,ct.Z)))
ax.set_xlabel('x')
ax.set_ylabel('z')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(ct.X, ct.Z, np.real(ct.b_par_ana(ct.X,ct.Z)))
ax.set_xlabel('x')
ax.set_ylabel('z')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(ct.X, ct.Z, np.imag(ct.b_par_ana(ct.X,ct.Z)))
ax.set_xlabel('x')
ax.set_ylabel('z')