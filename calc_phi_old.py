from scipy.integrate import solve_bvp
from scipy.io import FortranFile

def vA(z):
	return (g.epsilon + (1 - g.epsilon) * (1 + np.tanh(np.sin(np.pi * z / Lz) / g.h0)) / 2)

def wave_eqn(z, phi, p):
	omega = p[0]
	# d phi[0] / dt = phi[1]
	# d phi[1] / dt = -(omega / vA) ^ 2 phi[0]
	return np.vstack((phi[1], -(omega / vA(z) / np.cos(g.alpha)) ** 2 * phi[0]))

def bcs(phi_a, phi_b, p):
	# Impose periodic BCs, where phi(z=a) = phi(z=b) = 1
	# d phi / dz | z=a = d phi /dz | z=b
	omega = p[0]
	return np.array([phi_a[0]-1, phi_b[0]-1, phi_a[1] - phi_b[1]])

z_max = g.zb[-2]
z_min = g.zb[0]
Lz = (z_max - z_min) / 2

z = np.linspace(z_min, z_max, 9)
phi = np.zeros((2, z.size))
phi[0, 0] = 1
sol = solve_bvp(wave_eqn, bcs, z, phi, p=[2])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(g.zb[0:-1], vA(g.zb[0:-1]))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(g.zb[0:-1], sol.sol(g.zb[0:-1])[0])

f = FortranFile('Saves/omega.dat', 'w')
f.write_record(sol.p[0])
f.close()

f = FortranFile('Saves/phi.dat', 'w')
f.write_record(sol.sol(g.zb[1:-1])[0]) # phi_b
f.write_record(sol.sol(g.zb[1:-1])[1]) # dphi_b_dz
f.write_record(sol.sol(g.zc[1:-1])[0]) # phi_c
f.write_record(sol.sol(g.zc[1:-1])[1]) # dphi_c_dz
f.close()
