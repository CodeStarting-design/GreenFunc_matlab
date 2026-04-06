function G_A_xx = calculate_GAxx_SDP_FLAM(valid_poles,rho, h, er, freq) 
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);

    D = -1i * mu_0 / (8 * pi);

    Sum_R = 0;
    for p = 1:length(valid_poles)
        kp = valid_poles(p);
        u_p = h * sqrt(k1^2 - kp^2);
        
        if abs(u_p) < 1e-10 % 奇异性保护
            continue;
        end
        
        v_p = -u_p * cot(u_p);
        
        % Eq 12a-12c
        P_kp = 1i * 2 * sin(u_p) * kp * h;
        dQ_dk = (-kp * h^2 / u_p) * (cos(u_p) - u_p*sin(u_p) - (u_p/v_p)*sin(u_p) + v_p*cos(u_p));
        
        Res_p = besselh(0, 2, kp * rho) * (P_kp / dQ_dk);
        Sum_R = Sum_R + Res_p;
    end

    % =====================================================
    % 步骤 4: 沿最速下降路径计算积分 (SDP)
    % =====================================================
%     disp('正在计算最速下降路径 (SDP) 积分...');
    f_gamma1 = @(x) compute_integrand(x, rho, k0, k1, h, 1);
    f_gamma2 = @(x) compute_integrand(x, rho, k0, k1, h, 2);
    
    % 积分区间 [-inf, 0] 和 [0, inf]
    I_gamma1 = integral(f_gamma1, -Inf, 0, 'RelTol', 1e-4, 'AbsTol', 1e-6);
    I_gamma2 = integral(f_gamma2, 0, Inf, 'RelTol', 1e-4, 'AbsTol', 1e-6);
    % G_Gamma可以预计算（rho>0.05*lambda，后插值来提速）
    G_Gamma = D * (I_gamma1 + I_gamma2);

    % =====================================================
    % 步骤 5: 汇总结果 (Eq 11)
    % =====================================================

    G_A_xx = G_Gamma - 1i * 2 * pi * D * Sum_R;
end

% =========================================================
% 辅助函数 1: 积分核函数 (处理不同黎曼面)
% =========================================================
function val = compute_integrand(x, rho, k0, k1, h, sheet)
    % k_rho 变换 (Eq 14下方)
    k_rho = k0 - 1i * (x.^2) / (2 * rho);
    H02 = besselh(0, 2, k_rho * rho);
    
    % 区分黎曼面
    kz0_base = sqrt(k0^2 - k_rho.^2);
    if sheet == 1
        kz0 = kz0_base;       % 顶层黎曼面 (R1)
    else
        kz0 = -kz0_base;      % 底层黎曼面 (R2)
    end
    
    kz1 = sqrt(k1^2 - k_rho.^2);
    
    % 分母: j*k0z + k1z*cot(k1z*h)
    denom = 1i * kz0 + kz1 .* cot(kz1 * h);
    
    % Eq 14b, 14c 核心部分
    val = (x / rho) .* (2 * k_rho ./ denom) .* H02;
end