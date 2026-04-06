import numpy as np
from scipy.integrate import quad

def sommerfeld(rho, z, k0):
    def integrand(krho):
        kz0 = np.sqrt(k0**2 - krho**2 + 0j)
        if kz0.imag > 0:
            kz0 = kz0.real - 1j * kz0.imag
        from scipy.special import j0
        return (krho / kz0) * j0(krho * rho) * np.exp(-1j * kz0 * z)
    
    # complex integration is tricky, just check real part for small k0
    pass

