freq = 8e9; c = 3e8; lambda0 = c / freq; h = lambda0; er = 4;
mu_0 = 4 * pi * 1e-7; eps_0 = 8.854e-12; omega = 2 * pi * freq;
k0 = omega * sqrt(mu_0 * eps_0); k1 = omega * sqrt(mu_0 * er * eps_0);

valid_poles = FLAM_TE(0.1, h, er, freq);
rho = 0.1;

% From SDP_FLAM
Sum_R_SDP = 0;
for p = 1:length(valid_poles)
    kp = valid_poles(p);
    u_p = h * sqrt(k1^2 - kp^2);
    v_p = -u_p * cot(u_p);
    P_kp = 1i * 2 * sin(u_p) * kp * h;
    dQ_dk = (-kp * h^2 / u_p) * (cos(u_p) - u_p*sin(u_p) - (u_p/v_p)*sin(u_p) + v_p*cos(u_p));
    Res_p = besselh(0, 2, kp * rho) * (P_kp / dQ_dk);
    Sum_R_SDP = Sum_R_SDP + Res_p;
end

% From DCIM2
Sum_R_DCIM = 0;
for p = 1:length(valid_poles)
    kp = valid_poles(p);
    delta_k = 1e-8 * kp;
    
    % get_G_total
    krho = kp + delta_k;
    k0z = sqrt(k0^2 - krho^2);
    k1z = sqrt(k1^2 - krho^2);
    k0z(imag(k0z) > 0) = real(k0z(imag(k0z) > 0)) - 1j * abs(imag(k0z(imag(k0z) > 0)));
    k1z(imag(k1z) > 0) = real(k1z(imag(k1z) > 0)) - 1j * abs(imag(k1z(imag(k1z) > 0)));
    E = exp(-2 * 1i * k1z * h);
    Numerator = 2 * (1 - E);
    Denominator = 1i * ((k0z + k1z) - (k0z - k1z) * E);
    Gp = Numerator / Denominator;
    
    Res_G = delta_k * Gp;
    Res_p = kp * Res_G * besselh(0, 2, kp * rho);
    Sum_R_DCIM = Sum_R_DCIM + Res_p;
end

disp(['Sum_R_SDP = ', num2str(Sum_R_SDP)]);
disp(['Sum_R_DCIM = ', num2str(Sum_R_DCIM)]);
disp(['Ratio = ', num2str(Sum_R_DCIM / Sum_R_SDP)]);
