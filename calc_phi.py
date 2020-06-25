from scipy.integrate import solve_bvp
import control as ct

def wave_eqn(z, phi, p):
	omega = p[0]
	# d phi[0] / dt = phi[1]
	# d phi[1] / dt = -(omega / vA) ^ 2 phi[0]
	return np.vstack((phi[1], -(omega / vA(z) / np.cos(g.alpha)) ** 2 * phi[0]))