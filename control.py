import numpy as np

alpha = 0.25 * np.pi
k_perp = 2 * np.pi
omega_r = np.pi * np.cos(alpha)
omega_i = 0.1
omega = omega_r + 1j * omega_i
vA0 = 1

def vA(x,z):
	return vA0

def ux_ana(x,z):
	return ux0 * np.exp(1j * (kx * x + kz * z))

def b_par_ana(x,z):
	return b_par0 * np.exp(1j * (kx * x + kz * z))

xi = -2 * omega_i / omega_r

nx = 128
lx = 8 * abs(xi)
x_min = -lx
x_max =  lx
dx = (x_max - x_min) / (nx - 1)
x = np.linspace(x_min, x_max, nx)

nz = 128
Lz = 1
z_min = 0
z_max =  2 * Lz
dz = (z_max - z_min) / (nz - 1)
z = np.linspace(z_min, z_max, nz)

X, Z = np.meshgrid(x, z)

kx = np.pi / lx
kz = 2 * np.pi / Lz
ux0 = 1
b_par0 = (kz ** 2 * np.cos(alpha) ** 2 - (omega / vA0) ** 2) * ux0 / omega / kx

ux_x_min    = ux_ana(x_min, z)
b_par_x_min = b_par_ana(x_min, z)