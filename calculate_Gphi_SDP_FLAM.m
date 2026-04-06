function G_V = calculate_Gphi_SDP_FLAM(valid_poles, rho, h, er, freq) 
    % =====================================================
    % 常数与基本参数定义
    % =====================================================
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);

    % G_V 公式 18 的前置系数 D_v
    D_v = -1i / (8 * pi * eps_0);

    % =====================================================
    % 步骤 1: 计算极点留数之和 (Sum_R) 
    % 对应文献图 3 下方公式
    % =====================================================
    Sum_R = 0;
    for p = 1:length(valid_poles)
        kp = valid_poles(p);
        
        % 根据文献 Eq 4, 5 下方的定义计算 u 和 v
        u_p = h * sqrt(k1^2 - kp^2);
        
        if abs(u_p) < 1e-10 % 奇异性保护
            continue;
        end
        
        % 确定极点处 kz0 的分支。对表面波而言 kz0 是负虚数
        kz0_p = sqrt(k0^2 - kp^2);
        if imag(kp) >= 0
            kz0_p = -kz0_p;
        end
        v_p = 1i * kz0_p * h;
        
        % -----------------------------------------------------
        % 留数分子 N(kp) 
        % 对应文献中 N(k_rho_p) = (k0^2*v^2*er + k0z^2*u^2)*sin(2u) - 2*kp^2*v*u*sin^2(u)
        % -----------------------------------------------------
        N_kp = (k0^2 * v_p^2 * er + kz0_p^2 * u_p^2) * sin(2 * u_p) - ...
               2 * kp^2 * v_p * u_p * sin(u_p)^2;
        
        % -----------------------------------------------------
        % 留数分母的导数 dQ/dk (即 dD/dk_rho)
        % -----------------------------------------------------
        dQ_dk = 1i * (kp^2*h*v_p) * ( ...
                (v_p/u_p - u_p/v_p) * (er * cos(u_p)^2 - sin(u_p)^2) + ...
                cos(2 * u_p) * (er * v_p^2 / u_p - u_p) - ...
                sin(2 * u_p) * (er + 1) * (v_p + 1) ...
                );
        
        % 计算单个极点的留数并累加
        Res_p = besselh(0, 2, kp * rho) * (N_kp / dQ_dk);
        Sum_R = Sum_R + Res_p;
    end

    % =====================================================
    % 步骤 2: 沿最速下降路径计算连续谱积分 (SDP)
    % 对应公式 18 以及你提供的 Phi_V 截图
    % =====================================================
    f_gamma1 = @(x) compute_integrand_GV(x, rho, k0, k1, h, er, 1);
    f_gamma2 = @(x) compute_integrand_GV(x, rho, k0, k1, h, er, 2);
    
% x=-5:0.0001:5;
%     for i=1:length(x)
%         f1(i)=exp(x(i)*x(i)/2)*compute_integrand_GV(x(i), rho, k0, k1, h, er, 1);
%         f2(i)=exp(x(i)*x(i)/2)*compute_integrand_GV(x(i), rho, k0, k1, h, er, 2);
%     end
%     
%     figure
%     plot(x,f1)
%     figure
%     plot(x,f2)

    % 积分区间 [-inf, 0] 和 [0, inf]
    I_gamma1 = integral(f_gamma1, -Inf, 0, 'RelTol', 1e-4, 'AbsTol', 1e-6);
    I_gamma2 = integral(f_gamma2, 0, Inf, 'RelTol', 1e-4, 'AbsTol', 1e-6);
    
    G_Gamma = D_v * (I_gamma1 + I_gamma2);

    % =====================================================
    % 步骤 3: 汇总结果 (Eq 14 逻辑)
    % =====================================================
    G_V = G_Gamma - 1i * 2 * pi * D_v * Sum_R;
end

% =========================================================
% 辅助函数: G_V 积分核函数 
% =========================================================
function val = compute_integrand_GV(x, rho, k0, k1, h, er, sheet)
    % k_rho 变换
    k_rho = k0 - 1i * (x.^2) / (2 * rho);
    H02 = besselh(0, 2, k_rho * rho);
    
    % 区分黎曼面计算 kz0
    kz0_base = sqrt(k0^2 - k_rho.^2);
    if sheet == 1
        kz0 = kz0_base;       % 顶层黎曼面 (R1)
    else
        kz0 = -kz0_base;      % 底层黎曼面 (R2)
    end
    
    kz1 = sqrt(k1^2 - k_rho.^2);
    
    % 定义 u 和 v 
    u = kz1 * h;
    v = 1i * kz0 * h; 
    
    % 根据你提供的截图公式计算 Phi_V(x)
    % (exp(x^2/2) 与 公式18 积分表达式的 exp(-x^2/2) 抵消，不予体现)
    term_num = 2*(k0^2 .* v.^2 .* er + kz0.^2 .* u.^2) .* cot(u) - 2 .* v .* u .* k_rho.^2;
    term_den = kz0 .* k_rho .* (v + u .* cot(u)) .* (er .* v .* cot(u) - u);
    
    val = -1i .* (x ./ rho) .* (term_num ./ term_den) .* H02;
end