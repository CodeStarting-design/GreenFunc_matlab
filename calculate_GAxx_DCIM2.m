function G_A_xx = calculate_GAxx_DCIM2(valid_poles, rho, h, er, freq) 
    % =====================================================
    % 常数与基本参数初始化
    % =====================================================
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);

    % 全局系数 D = -j*mu0/(8*pi)
    D = -1i * mu_0 / (8 * pi);
    valid_poles = valid_poles(:); 
    
    % =====================================================
    % 步骤 1: 数值留数提取 (对 get_G_total 的留数)
    % =====================================================
    Sum_R = 0;
    Res_G = zeros(length(valid_poles), 1);
    
    for p = 1:length(valid_poles)
        kp = valid_poles(p);
        
        % Res_G = lim_{k->kp} (k - kp) * get_G_total(k)
        delta_k = 1e-6 * kp;
        Gp = get_G_total(kp + delta_k, k0, k1, h);
        Res_G(p) = delta_k * Gp;
        
        % 空间域极点贡献: Res[F * H0^(2)] = kp * Res_G * H0^(2)
        Res_p = kp * Res_G(p) * besselh(0, 2, kp * rho);
        Sum_R = Sum_R + Res_p;
    end

    % =====================================================
    % 步骤 2: 定义 DCIM 采样路径 (沿 kz0 虚轴)
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
    % 步骤 3: 构建谱域尾项 (准静态+极点减除+GPOF拟合)
    % =====================================================
    G_tot_path = get_G_total(krho_path, k0, k1, h);
    G_qs_path = (1 - exp(-2 * 1i * kz0_path * h)) ./ (1i .* kz0_path);
    
    % 极点部分分式: get_G_total 的极点展开 = 2*kp*Res_G/(krho^2-kp^2)
    G_swp_path = zeros(size(krho_path));
    for i = 1:length(krho_path)
        diff_poles = krho_path(i).^2 - valid_poles.^2;
        diff_poles(abs(diff_poles) < 1e-25) = 1e-25; 
        G_swp_path(i) = sum((2 .* valid_poles .* Res_G) ./ diff_poles);
    end
    
    Ftail = G_tot_path - G_qs_path - G_swp_path;
    
    % 乘以 j*kz0 使 Ftail 有界，方便 GPOF 拟合
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
    
    % 映射回物理空间: 复数图像深度 z_i = -alpha_t/k0
    a_DCIM = a_t .* exp(-T02 .* alpha_t);
    alpha_DCIM = -1j * alpha_t ./ k0;

    % =====================================================
    % 步骤 5: 索末菲恒等式空间域转换
    % =====================================================
    % 索末菲恒等式: int[exp(-jkz0*z)/kz0] * krho * H02 dkrho = -2i*exp(-jk0R)/R
    % 对 a/(jkz0) 形式: 系数 = (-2i)/j = +2, 含 SDP 雅可比 i 后 = +2i
    
    % 1. GPOF 图像项
    I_DCIM = 0;
    for i = 1:length(a_DCIM)
        Rc = sqrt(rho^2 - alpha_DCIM(i)^2); 
        I_DCIM = I_DCIM + (2i) * a_DCIM(i) * (exp(-1j * k0 * Rc) / Rc);
    end

    % 2. 准静态项: G_qs = (1-exp(-2jkz0*h))/(jkz0) -> 两个球面波
    R0 = rho;
    R1 = sqrt(rho.^2 + (2*h).^2);
    I_qs = (2i) * (exp(-1j * k0 * R0) / R0 - exp(-1j * k0 * R1) / R1);

    % =====================================================
    % 步骤 6: 汇总 (对标 SDP-FLAM: G = D*(I_gamma - 2*pi*i*Sum_R))
    % =====================================================
    I_total = I_DCIM + I_qs - 2i * pi * Sum_R;
    G_A_xx = D * I_total;

end

% =========================================================
% 谱域函数: get_G_total = 2*K11 = 2*sin(kz1*h)/[kz1*cos(kz1*h)+i*kz0*sin(kz1*h)]
% =========================================================
function G = get_G_total(krho, k0, k1, h)
    k0z = sqrt(k0.^2 - krho.^2);
    k1z = sqrt(k1.^2 - krho.^2);
    
    % 强制辐射/衰减条件 (Im(kz) < 0)
    k0z(imag(k0z) > 0) = real(k0z(imag(k0z) > 0)) - 1j * abs(imag(k0z(imag(k0z) > 0)));
    k1z(imag(k1z) > 0) = real(k1z(imag(k1z) > 0)) - 1j * abs(imag(k1z(imag(k1z) > 0)));

    E = exp(-2 * 1i * k1z * h);
    Numerator = 2 * (1 - E);
    Denominator = 1i * ((k0z + k1z) - (k0z - k1z) .* E);
    G = Numerator ./ Denominator;
end

% =========================================================
% GPOF 独立模块 (广义铅笔函数法)
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
