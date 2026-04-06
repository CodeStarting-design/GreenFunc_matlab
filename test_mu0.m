f = 8e9; c = 3e8; lambda0 = c / f; h = lambda0; er = 4;
valid_poles1 = FLAM_TE(0.1, h, er, f);
disp('SDP-FLAM:');
calculate_GAxx_SDP_FLAM(valid_poles1, 0.1, h, er, f)
disp('DCIM2:');
calculate_GAxx_DCIM2(valid_poles1, 0.1, h, er, f)
