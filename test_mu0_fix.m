f = 8e9; c = 3e8; lambda0 = c / f; h = lambda0; er = 4;
mu_0 = 4 * pi * 1e-7;
valid_poles1 = FLAM_TE(0.1, h, er, f);

G_SDP = calculate_GAxx_SDP_FLAM(valid_poles1, 0.1, h, er, f);

% DCIM without mu_0
G_DCIM = calculate_GAxx_DCIM2(valid_poles1, 0.1, h, er, f) * mu_0; 
% wait, DCIM divides by mu_0, so multiply back

disp(['G_SDP = ', num2str(G_SDP)]);
disp(['G_DCIM = ', num2str(G_DCIM)]);
disp(['Direct = 209.704913930695 - 20.8135378016034i']); % rho=0.1 might not be this direct value, let's just see G_SDP vs G_DCIM
