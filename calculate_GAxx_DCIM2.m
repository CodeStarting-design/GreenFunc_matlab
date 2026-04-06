function G_A_xx = calculate_GAxx_DCIM2(valid_poles, rho, h, er, freq) 
    % =====================================================
    % 常数与基本参数初始化
    % =====================================================
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);

    D = -1i * mu_0 / (8 * pi);
    valid_poles = valid_poles(:); 
    
    % =====================================================
    % 步骤 1: 数值留数
    % =====================================================
    Sum_R = 0;
    Res_G = zeros(length(valid_poles), 1);
    
    for p = 1:length(valid_poles)
        kp = valid_poles(p);
        
        % Res_G = lim_{k->kp} (k-kp) * get_G_total(k)
        delta_k = 1e-8 * kp;
        Gp = get_G_total(kp + delta_k, k0, k1, h);
        Res_G(p) = delta_k * Gp;
        
        % 空间域极点: Res[F*H02] at kp = kp*Res_G*H02
        Res_p = kp * Res_G(p) * besselh(0, 2, kp * rho);
        Sum_R = Sum_R + Res_p;
    end

    % =====================================================
    % 步骤 2: DCIM 采样路径
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
    % 步骤 3: 谱域尾项 (用数值留数做极点减除)
    % =====================================================
    G_tot_path = get_G_total(krho_path, k0, k1, h);
    G_qs_path = (1 - exp(-2 * 1i * kz0_path * h)) ./ (1i .* kz0_path);
    
    G_swp_path = zeros(size(krho_path));
    for i = 1:length(krho_path)
        diff_poles = krho_path(i).^2 - valid_poles.^2;
        diff_poles(abs(diff_poles) < 1e-25) = 1e-25; 
        G_swp_path(i) = sum((2 .* valid_poles .* Res_G) ./ diff_poles);
    end
    
    Ftail = G_tot_path - G_qs_path - G_swp_path;
    y_fit = Ftail .* (1j * kz0_path);
    
    if any(isnan(y_fit)) || any(isinf(y_fit))
        error('DCIM 采样路径出现 NaN 或 Inf');
    end

    % =====================================================
    % 步骤 4: GPOF 拟合
    % =====================================================
    tol_svd = 1e-4; 
    tol_eig = 1e-16; 
    max_images = 10;
    
    [a_t, alpha_t] = RunGPOF_Standalone(y_fit, dt, tol_svd, tol_eig, max_images);
    
    a_DCIM = a_t .* exp(-T02 .* alpha_t);
    alpha_DCIM = -1j * alpha_t ./ k0;

    % =====================================================
    % 步骤 5: 空间域 (索末菲恒等式)
    % =====================================================
    % 索末菲恒等式 (Chew): int[exp(-jkz0*z)/kz0]*krho*H02 dkrho = +2i*exp(-jk0R)/R
    % 对 a/(j*kz0) 形式: 系数 = +2i/j = +2
    % I_gamma = i * SIP积分 => 含雅可比后系数 = i*(+2) = +2i

    % GPOF 图像项
    I_DCIM = 0;
    for i = 1:length(a_DCIM)
        Rc = sqrt(rho^2 - alpha_DCIM(i)^2); 
        I_DCIM = I_DCIM + (2i) * a_DCIM(i) * (exp(-1j * k0 * Rc) / Rc);
    end

    % 准静态项
    R0 = rho;
    R1 = sqrt(rho.^2 + (2*h).^2);
    I_qs = (2i) * (exp(-1j * k0 * R0) / R0 - exp(-1j * k0 * R1) / R1);

    % =====================================================
    % 步骤 6: 汇总 (与 SDP-FLAM 对标)
    % SDP: G = D*(I_gamma - 2*pi*i*Sum_R)
    % DCIM: I_cont 替代 I_gamma, Sum_R 用相同的解析留数
    % =====================================================
    I_total = I_DCIM + I_qs + 2 * pi * Sum_R;
    G_A_xx = D * I_total / mu_0;

end

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
