import numpy as np
from scipy.special import hankel2
from scipy.optimize import fsolve

freq = 8e9
c = 3e8
h = c / freq
er = 4
mu_0 = 4 * np.pi * 1e-7
eps_0 = 8.854e-12
omega = 2 * np.pi * freq
k0 = omega * np.sqrt(mu_0 * eps_0)
k1 = omega * np.sqrt(mu_0 * er * eps_0)

def denom_func(krho):
    k0z = np.sqrt(k0**2 - krho**2 + 0j)
    if k0z.imag > 0:
        k0z = k0z.real - 1j * abs(k0z.imag)
    k1z = np.sqrt(k1**2 - krho**2 + 0j)
    if k1z.imag > 0:
        k1z = k1z.real - 1j * abs(k1z.imag)
    
    E = np.exp(-2j * k1z * h)
    denom = 1j * ((k0z + k1z) - (k0z - k1z) * E)
    return denom.real, denom.imag

# k0 = 167.6, k1 = 335.2
# Let's search for roots in (k0, k1)
k_test = np.linspace(k0+1, k1-1, 1000)
denoms = [denom_func(k) for k in k_test]
abs_denoms = [np.abs(d[0] + 1j*d[1]) for d in denoms]

# find local minima
poles = []
for i in range(1, len(abs_denoms)-1):
    if abs_denoms[i] < abs_denoms[i-1] and abs_denoms[i] < abs_denoms[i+1]:
        from scipy.optimize import root_scalar
        res = root_scalar(lambda k: denom_func(k)[1], bracket=[k_test[i-1], k_test[i+1]])
        if abs(denom_func(res.root)[0]) < 1e-5:
            poles.append(res.root)

print("Found poles:", poles)

for kp in poles:
    rho = 0.1
    # DCIM2 formulation
    delta_k = 1e-8 * kp
    krho = kp + delta_k
    k0z = np.sqrt(k0**2 - krho**2 + 0j)
    if k0z.imag > 0: k0z = k0z.real - 1j * abs(k0z.imag)
    k1z = np.sqrt(k1**2 - krho**2 + 0j)
    if k1z.imag > 0: k1z = k1z.real - 1j * abs(k1z.imag)
    E = np.exp(-2j * k1z * h)
    Numerator = 2 * (1 - E)
    Denominator = 1j * ((k0z + k1z) - (k0z - k1z) * E)
    Gp = Numerator / Denominator
    Res_G = delta_k * Gp
    Res_p_DCIM = kp * Res_G * hankel2(0, kp * rho)
    
    # SDP formulation
    u_p = h * np.sqrt(k1**2 - kp**2 + 0j)
    v_p = -u_p / np.tan(u_p)
    P_kp = 1j * 2 * np.sin(u_p) * kp * h
    dQ_dk = (-kp * h**2 / u_p) * (np.cos(u_p) - u_p*np.sin(u_p) - (u_p/v_p)*np.sin(u_p) + v_p*np.cos(u_p))
    Res_p_SDP = hankel2(0, kp * rho) * (P_kp / dQ_dk)
    
    print(f"Pole: {kp}")
    print(f"Res_p_DCIM: {Res_p_DCIM}")
    print(f"Res_p_SDP: {Res_p_SDP}")
    print(f"Ratio: {Res_p_DCIM / Res_p_SDP}")

