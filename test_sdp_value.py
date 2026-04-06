import numpy as np
from scipy.integrate import quad
from scipy.special import hankel2

freq = 8e9
c = 3e8
h = c / freq
er = 4
mu_0 = 4 * np.pi * 1e-7
eps_0 = 8.854e-12
omega = 2 * np.pi * freq
k0 = omega * np.sqrt(mu_0 * eps_0)
k1 = omega * np.sqrt(mu_0 * er * eps_0)
D = -1j * mu_0 / (8 * np.pi)

rho = 0.01 * h

def compute_integrand(x, sheet):
    k_rho = k0 - 1j * (x**2) / (2 * rho)
    H02 = hankel2(0, k_rho * rho)
    kz0_base = np.sqrt(k0**2 - k_rho**2 + 0j)
    if sheet == 1:
        kz0 = kz0_base
    else:
        kz0 = -kz0_base
    kz1 = np.sqrt(k1**2 - k_rho**2 + 0j)
    denom = 1j * kz0 + kz1 / np.tan(kz1 * h)
    val = (x / rho) * (2 * k_rho / denom) * H02
    return val

def f_real(x): return compute_integrand(x, 1).real
def f_imag(x): return compute_integrand(x, 1).imag

# Just evaluate at x=0.1 to see the scale
v = compute_integrand(0.1, 1)
print("Integrand at x=0.1:", v)
print("D scale:", D)
print("D * Integrand:", D * v)

