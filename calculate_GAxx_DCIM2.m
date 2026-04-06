function G_A_xx = calculate_GAxx_DCIM2(valid_poles, rho, h, er, freq) 
    % =====================================================
    % 常数与基本参数初始化
    % =====================================================
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);

    % 全局系数 D
    D = -1i * mu_0 / (8 * pi);
    valid_poles = valid_poles(:); 
    
    % =====================================================
    % 步骤 1: 绝对精确的数值留数提取 (彻底抛弃解析求导的风险)
    % =====================================================
    Sum_R = 0;
    Res_G = zeros(length(valid_poles), 1); % get_G_total 的数值留数
    
    for p = 1:length(valid_poles)
        kp = valid_poles(p);
        
        % 数值留数: Res_G = lim_{k->kp} (k - kp) * get_G_total(k)
        delta_k = 1e-6 * kp;
        Gp = get_G_total(kp + delta_k, k0, k1, h);
        Res_G(p) = delta_k * Gp;
        
        % 同时计算解析留数做对比
        u_p = h * sqrt(k1^2 - kp^2);
        if abs(u_p) > 1e-10
            v_p = -u_p * cot(u_p);
            P_kp = 1i * 2 * sin(u_p) * kp * h;
            dQ_dk = (-kp * h^2 / u_p) * (cos(u_p) - u_p*sin(u_p) - (u_p/v_p)*sin(u_p) + v_p*cos(u_p));
            Res_i_analytic = P_kp / dQ_dk;
            fprintf('  pole kp/k0=%.4f: Res_G=%.4e%+.4ej, kp*Res_G=%.4e%+.4ej, Res_i=%.4e%+.4ej, ratio=%.4f%+.4fj\n', ...
                real(kp/k0), real(Res_G(p)), imag(Res_G(p)), ...
                real(kp*Res_G(p)), imag(kp*Res_G(p)), ...
                real(Res_i_analytic), imag(Res_i_analytic), ...
                real(Res_i_analytic/(kp*Res_G(p))), imag(Res_i_analytic/(kp*Res_G(p))));
        end
        
        % 空间域极点贡献: F(kp)的留数 = kp * Res_G
        Res_p = kp * Res_G(p) * besselh(0, 2, kp * rho);
        Sum_R = Sum_R + Res_p;
    end
   

    % =====================================================
    % 步骤 2: 定义 DCIM 采样路径
    % =====================================================
    N_samples = 300;                  
    T02 = 1.2;          
    T01 = 600.0;                     
    t = linspace(0.0, T01, N_samples).'; 
    dt = t(2) - t(1);
    
    kz0_path = -1j * k0 * (T02 + t);
    krho_path = sqrt(k0^2 - kz0_path.^2);
    krho_path = abs(real(krho_path)) + 1j * abs(imag(krho_path));

    % =====================================================
    % 步骤 3: 构建谱域尾项 (平滑化)
    % =====================================================
    G_tot_path = get_G_total(krho_path, k0, k1, h);
    G_qs_path = (1 - exp(-2 * 1i * kz0_path * h)) ./ (1i .* kz0_path);
    
    % 减去匹配的极点
    G_swp_path = zeros(size(krho_path));
    for i = 1:length(krho_path)
        diff_poles = krho_path(i).^2 - valid_poles.^2;
        diff_poles(abs(diff_poles) < 1e-25) = 1e-25; 
        G_swp_path(i) = sum((2 .* valid_poles .* Res_G) ./ diff_poles);
    end
    
    Ftail = G_tot_path - G_qs_path - G_swp_path;
    
    % 拟合目标: Ftail * (j*kz0) 
    y_fit = Ftail .* (1j * kz0_path);
    
    if any(isnan(y_fit)) || any(isinf(y_fit))
        error('DCIM 采样路径出现 NaN 或 Inf，请检查参数！');
    end

    % =====================================================
    % 步骤 4: 执行 GPOF 拟合
    % =====================================================
    tol_svd = 1e-4; 
    tol_eig = 1e-16; 
    max_images = 10;
    
    [a_t, alpha_t] = RunGPOF_Standalone(y_fit, dt, tol_svd, tol_eig, max_images);
    
    % 映射回 kz0 空间
    a_DCIM = a_t .* exp(-T02 .* alpha_t);
    alpha_DCIM = -1j * alpha_t ./ k0;

    % =====================================================
    % 步骤 5: 利用索末菲恒等式计算空间域积分
    % =====================================================
    
    % --- 诊断: 直接数值积分 Ftail 在 SIP 上的贡献 ---
    % 对 Ftail * krho * H0^(2) 做数值积分 (梯形法则)
    % 这绕过 GPOF，直接验证 Ftail 谱域分解的正确性
    dkrho = diff(krho_path);
    integrand_tail = Ftail .* krho_path .* besselh(0, 2, krho_path * rho);
    I_tail_direct = sum(0.5*(integrand_tail(1:end-1)+integrand_tail(2:end)) .* dkrho);
    
    % 准静态部分同样直接积分
    integrand_qs = G_qs_path .* krho_path .* besselh(0, 2, krho_path * rho);
    I_qs_direct = sum(0.5*(integrand_qs(1:end-1)+integrand_qs(2:end)) .* dkrho);
    
    % 极点部分同样直接积分
    integrand_swp = G_swp_path .* krho_path .* besselh(0, 2, krho_path * rho);
    I_swp_direct = sum(0.5*(integrand_swp(1:end-1)+integrand_swp(2:end)) .* dkrho);
    
    % 连续谱总积分 = tail + qs + swp 在 SIP 上的贡献
    I_cont_SIP = I_tail_direct + I_qs_direct + I_swp_direct;
    
    % I_gamma = i * I_SIP (SDP 雅可比)
    % 但 DCIM 路径不是 SDP 路径！DCIM 路径是沿虚轴的采样路径。
    % 尝试两种方式，看哪个匹配 SDP:
    
    % 方式A: 直接把 SIP 积分结果视为 I_gamma 的一部分
    % G = (1/(8*pi)) * [I_cont_SIP + pole_contrib]
    % 其中 pole_contrib = -2*pi*i * Sum_R
    G_direct_8pi = (1/(8*pi)) * (I_cont_SIP - 2i*pi*Sum_R);
    
    % 方式B: 用 D 系数
    G_direct_D = D * (I_cont_SIP - 2i*pi*Sum_R);
    
    % 方式C: 用 D*i 系数 (SDP 雅可比)
    G_direct_Di = D * 1i * (I_cont_SIP - 2i*pi*Sum_R);
    
    % 打印诊断信息
    fprintf('rho=%.4e: I_SIP=%.6e%+.6ej, Sum_R=%.6e%+.6ej\n', ...
        rho, real(I_cont_SIP), imag(I_cont_SIP), real(Sum_R), imag(Sum_R));
    fprintf('  G_1/(8pi) = %.6e%+.6ej\n', real(G_direct_8pi), imag(G_direct_8pi));
    fprintf('  G_D       = %.6e%+.6ej\n', real(G_direct_D), imag(G_direct_D));
    fprintf('  G_D*i     = %.6e%+.6ej\n', real(G_direct_Di), imag(G_direct_Di));
    
    % 1. GPOF 空域积分
    I_DCIM = 0;
    for i = 1:length(a_DCIM)
        Rc = sqrt(rho^2 - alpha_DCIM(i)^2); 
        I_DCIM = I_DCIM + (-2i) * a_DCIM(i) * (exp(-1j * k0 * Rc) / Rc);
    end

    % 2. 准静态项空域积分
    R0 = rho;
    R1 = sqrt(rho.^2 + (2*h).^2);
    I_qs = (-2i) * (exp(-1j * k0 * R0) / R0 - exp(-1j * k0 * R1) / R1);

    % =====================================================
    % 步骤 6: 汇总最终结果
    % =====================================================
    I_total = I_DCIM + I_qs - 2i * pi * Sum_R;
    G_A_xx = D * I_total;
    
    fprintf('  G_DCIM(GPOF) = %.6e%+.6ej\n', real(G_A_xx), imag(G_A_xx));

end

% =========================================================
% 谱域函数
% =========================================================
function G = get_G_total(krho, k0, k1, h)
    k0z = sqrt(k0.^2 - krho.^2);
    k1z = sqrt(k1.^2 - krho.^2);
    
    k0z(imag(k0z) > 0) = real(k0z(imag(k0z) > 0)) - 1j * abs(imag(k0z(imag(k0z) > 0)));
    k1z(imag(k1z) > 0) = real(k1z(imag(k1z) > 0)) - 1j * abs(imag(k1z(imag(k1z) > 0)));

    E = exp(-2 * 1i * k1z * h);
    Numerator = 2 * (1 - E);
    Denominator = 1i * ((k0z + k1z) - (k0z - k1z) .* E);
    G = Numerator ./ Denominator;
end

% =========================================================
% GPOF 独立模块
% =========================================================
function [a, alpha] = RunGPOF_Standalone(y, dt, tol_svd, tol_eig, max_num_images)
    y = y(:); 
    N = length(y);
    L = floor(N / 2);
    if max_num_images > 0 && max_num_images < L
        L = max_num_images;
    end
    
    Y1 = zeros(N - L, L);
    Y2 = zeros(N - L, L);
    for ii = 1:L
        Y1(:, ii) = y(ii : ii + N - L - 1);
        Y2(:, ii) = y(ii + 1 : ii + N - L);
    end
    
    [U, S, V] = svd(Y1, 'econ');
    s_vals = diag(S);
    
    nst = find(s_vals < tol_svd * s_vals(1), 1) - 1;
    if isempty(nst), nst = length(s_vals); end
    if nst < 1, nst = 1; end
    
    U = U(:, 1:nst);
    S = S(1:nst, 1:nst);
    V = V(:, 1:nst);
    
    D_inv = diag(1 ./ diag(S));
    Z = D_inv * U' * Y2 * V;
    
    w = eig(Z);
    [~, sort_idx] = sort(abs(w), 'descend');
    w = w(sort_idx);
    
    nwt = find(abs(w) < tol_eig * abs(w(1)), 1) - 1;
    if isempty(nwt), nwt = length(w); end
    if nwt < 1, nwt = 1; end
    w = w(1:nwt);
    
    Y3 = zeros(N, nwt);
    for ii = 1:nwt
        Y3(:, ii) = w(ii).^(0 : N-1).';
    end
    
    a = Y3 \ y;
    alpha = log(w) / dt;
end