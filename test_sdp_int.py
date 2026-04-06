import numpy as np
from scipy.integrate import quad
from test_sdp_value import f_real, f_imag

res_r, _ = quad(f_real, 0, 100, limit=100)
res_i, _ = quad(f_imag, 0, 100, limit=100)

val = res_r + 1j * res_i
# multiply by 2 for -inf to inf
val *= 2

D = -5e-8j
print("I_gamma:", val)
print("D * I_gamma:", D * val)
print("D/mu_0 * I_gamma:", (-1j / (8*np.pi)) * val)
