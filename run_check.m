f = 8e9; c = 3e8; lambda0 = c / f; h = lambda0; er = 4;
mu_0 = 4 * pi * 1e-7; eps_0 = 8.854e-12;
valid_poles1 = FLAM_TE(0.1, h, er, f);

rho = 0.000375;
G_DCIM = calculate_GAxx_DCIM2(valid_poles1, rho, h, er, f) * mu_0;
disp(['G_DCIM at 0.01 lambda = ', num2str(G_DCIM)]);
