import numpy as np
import control as ct
import tridiag as td

def dUdx(U, t):
	ux = U[0:ct.nz-1]
	b_par = U[ct.nz:]
	u_perp = td.calc_u_perp(b_par, x)

	dux_dx[1:-1] = -(1j * ct.omega * ct.b_par[1:-1] + 1j * ct.k_perp * u_perp[1:-1] - \
					 np.sin(alpha) / (2 * ct.dz))
