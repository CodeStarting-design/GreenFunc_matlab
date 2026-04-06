import numpy as np
from scipy.special import hankel2
import cmath

freq = 8e9
c = 3e8
h = c / freq
er = 4
mu_0 = 4 * np.pi * 1e-7
eps_0 = 8.854e-12
omega = 2 * np.pi * freq
k0 = omega * np.sqrt(mu_0 * eps_0)
k1 = omega * np.sqrt(mu_0 * er * eps_0)

valid_poles = [224.283186358364] # example pole for TE mode at 8GHz, er=4, h=lambda0?
# Wait, let me extract the pole from FLAM_TE in python or just hardcode one.
# Let's extract the poles using the matlab output if possible.
