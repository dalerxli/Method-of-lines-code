import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import control as ct
import tridiag as td

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# surf = ax.plot_surface(ct.X, ct.Z, np.real(ct.ux_ana(ct.X,ct.Z)))
# ax.set_xlabel('x')
# ax.set_ylabel('z')

u_perp_x_min = td.calc_u_perp(ct.x_min,ct.b_par_x_min)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(ct.z, np.real(u_perp_x_min))
# ax.plot(ct.z, np.real(ct.u_perp_ana(ct.x_min,ct.z)))

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(ct.z, np.imag(u_perp_x_min))
# ax.plot(ct.z, np.imag(ct.u_perp_ana(ct.x_min,ct.z)))

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(ct.z, np.real(td.vector_d(ct.b_par_x_min)))
# ax.plot(ct.z, np.imag(td.vector_d(ct.b_par_x_min)))
