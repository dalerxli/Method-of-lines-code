import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
import control as ct

def matrix_A(x):
	a = np.full(ct.nz, np.cos(ct.alpha) ** 2 / ct.dz ** 2) # Off diagonal elements
	b = ct.omega ** 2 / ct.vA(x,ct.z) ** 2 - 2 * a         # Diagonal elements
	b[0]  += b[1]
	b[-1] += a[0] ** 2 / b[1]
	ab = np.zeros((3, ct.nz), dtype=complex)
	ab[0,:] = a
	ab[1,:] = b
	ab[2,:] = a
	return ab

def vector_d(b_par):
	d = np.zeros(ct.nz, dtype=complex)
	d[1:-1] = -ct.omega * (ct.k_perp * b_par[1:-1] + \
						   1j * np.sin(ct.alpha) / ct.dz * (b_par[2:] - b_par[0:-2]))
	d[0]    = -ct.omega * (ct.k_perp * b_par[0] + \
						   1j * np.sin(ct.alpha) / ct.dz * (b_par[1]  - b_par[-1]))
	d[-1]   = -ct.omega * (ct.k_perp * b_par[-1] + \
						   1j * np.sin(ct.alpha) / ct.dz * (b_par[0]  - b_par[-2]))
	return d

def calc_u_perp(x, b_par):

	ab = matrix_A(x)
	d = vector_d(b_par)

	a = ab[0,0]
	b1 = ab[1,1]
	u = np.zeros(ct.nz, dtype=complex)
	v = np.zeros(ct.nz, dtype=complex)
	u[0] = -b1
	u[-1] = a
	v[0] = 1
	v[-1] = -a / b1

	y = solve_banded((1,1), ab, d)
	q = solve_banded((1,1), ab, u)

	fig = plt.figure()
	plt.plot(np.real(y))
	plt.plot(np.imag(y))

	fig = plt.figure()
	plt.plot(np.real(q))
	plt.plot(np.imag(q))

	vTy = v[0] * y[0] + v[-1] * y[-1]
	vTq = v[0] * q[0] + v[-1] * q[-1]
	u_perp = y - vTy / (1 + vTq) * q

	return u_perp