function valid_poles = FLAM_TM(rho, h, er, freq)
    % --- 物理常数与波矢 ---
%     c = 3e8; 
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);
    
    % =====================================================
    % 步骤 1: FLAM - 寻找所有数学极点
    % =====================================================
%     disp('正在通过 FLAM 泰勒展开寻找 TM 模式复平面极点...');
    r2 = (k1^2 - k0^2) * h^2;
    N_max = 80; 
    a_points = -100:1:100; % 扫描实轴正半轴
    all_u_roots = [];
    
    for a = a_points
        % ---------------------------------------------------------
        % 构建 TM 模式多项式 (依据 Eq 2.102)
        % ---------------------------------------------------------
        
        % 1. 构造 Taylor 展开部分的多项式 P_taylor (最高次幂在前)
        P_taylor = zeros(1, N_max + 1);
        for n = 0:N_max
            P_taylor(N_max + 1 - n) = cos(a + n * pi / 2) / factorial(n);
        end
        
        % 2. 构造相乘的二次多项式 P_mult = [t^2系数, t^1系数, t^0系数]
        C2 = (1 - er^2) / 4;
        C1 = a * (1 - er^2) / 2;
        C0 = (1 - er^2) * a^2 / 4 + er^2 * r2;
        P_mult = [C2, C1, C0];
        
        % 两者相乘 (在多项式运算中即为卷积)
        P_part1 = conv(P_mult, P_taylor); % 结果长度为 N_max + 3
        
        % 3. 构造尾部的独立二次多项式 P_tail = [t^2系数, t^1系数, t^0系数]
        D2 = -(1 + er^2) / 4;
        D1 = -(1 + er^2) * a / 2;
        D0 = -(1 + er^2) * a^2 / 4 + er^2 * r2;
        P_tail = [D2, D1, D0];
        
        % 将 P_tail 补齐前导 0，使其长度与 P_part1 一致，以便相加
        P_tail_padded = [zeros(1, length(P_part1) - 3), P_tail];
        
        % 4. 合成最终的多项式系数
        coeffs = P_part1 + P_tail_padded;
        
        % 求解多项式根
        t_roots = roots(coeffs);
        u_roots = (t_roots + a) / 2;
        
        % ---------------------------------------------------------
        % 代回 TM 模式超越方程平方形式验证 (依据 Eq 2.100)
        % ---------------------------------------------------------
        for i = 1:length(u_roots)
            u_t = u_roots(i);
            % Eq 2.100: [(1-er^2)u^2 + er^2*r^2]cos(2u) - (er^2+1)u^2 + er^2*r^2 = 0
            val = ((1 - er^2)*u_t^2 + er^2*r2) * cos(2*u_t) - (er^2 + 1)*u_t^2 + er^2*r2;
            % 考虑到 TM 公式由于 er 的乘积项数值会被放大，容差适当放宽到 1e-7
            if abs(val) < 1e-7
                all_u_roots = [all_u_roots; u_t];
            end
        end
    end
    
    % --- 复数去重 ---
    tolerance = 1e-4;
    unique_u = [];
    for i = 1:length(all_u_roots)
        if isempty(unique_u)
            unique_u = [unique_u; all_u_roots(i)];
        else
            if min(abs(unique_u - all_u_roots(i))) > tolerance
                unique_u = [unique_u; all_u_roots(i)];
            end
        end
    end
    
    % 映射回 k_rho 平面
    k_rho_poles = sqrt(k1^2 - (unique_u ./ h).^2);

    % =====================================================
    % === k_rho 平面去重 ===
    % 消除 +u 和 -u 映射后产生的重复 k_rho
    % =====================================================
    unique_kp = [];
    for i = 1:length(k_rho_poles)
        if isempty(unique_kp)
            unique_kp = [unique_kp; k_rho_poles(i)];
        else
            % 计算距离，若大于容差则认为是新极点
            if min(abs(unique_kp - k_rho_poles(i))) > 1e-5
                unique_kp = [unique_kp; k_rho_poles(i)];
            end
        end
    end
    k_rho_poles = unique_kp; % 更新去重后的 k_rho

    % =====================================================
    % 步骤 2: 极点分类与拓扑剔除 (剔除增根)
    % =====================================================
    SWP_list = [];
    LWP_list = [];
    
    for i = 1:length(k_rho_poles)
        kp = k_rho_poles(i);
        
        % 基础物理范围限制 (向外传播，衰减)
        if real(kp) > 0 && real(kp) <= real(k1) * 1.1 && imag(kp) <= 1e-5
            
            % 1. 黎曼面 kz0 计算 (Top: 衰减, 强制负虚部)
            kz0_top = sqrt(k0^2 - kp^2);
            if imag(kz0_top) > 0, kz0_top = -kz0_top; end
            kz0_bot = -kz0_top; 
            
            kz1_val = sqrt(k1^2 - kp^2);
            
            % ---------------------------------------------------------
            % 2. 代回 TM 模式原始方程验算 (依据 Eq 2.98)
            % f(k_rho) = er * j * kz0 * h * cos(kz1*h) - kz1 * h * sin(kz1*h) = 0
            % 注意：为了避免代码出现 0，公式两边约掉公共乘子 h (因为 h 不等于 0)
            % ---------------------------------------------------------
            eq_top = er * 1i * kz0_top * cos(kz1_val * h) - kz1_val * sin(kz1_val * h);
            eq_bot = er * 1i * kz0_bot * cos(kz1_val * h) - kz1_val * sin(kz1_val * h);
            
            % 3. 分类判定
            if real(kp) > real(k0) && real(kp) < real(k1) && abs(imag(kp)) < 1e-4
                if abs(eq_top) < 1e-3  % 是真正的 TM 表面波极点
                    SWP_list = [SWP_list; kp];
                end
            else
                if abs(eq_bot) < 1e-3  % 是真正的 TM 漏波极点
                    LWP_list = [LWP_list; kp];
                end
            end
        end
    end
    
    fprintf('\n--- TM 模式物理极点识别结果 ---\n');
    fprintf('TM 表面波极点 (SWP) 共 %d 个:\n', length(SWP_list));
    for i = 1:length(SWP_list)
        fprintf('  SWP %d: kp/k0 = %.4f + %.4fi\n', i, real(SWP_list(i))/real(k0), imag(SWP_list(i))/real(k0));
    end
    fprintf('TM 漏波极点 (LWP) 共 %d 个\n', length(LWP_list));

    % 汇总物理极点
    physical_poles = [SWP_list; LWP_list];
    
    % --- 距离衰减判据 (洗刷远处的 LWP) ---
    valid_poles = physical_poles(abs(imag(physical_poles) .* rho) < 15);
end