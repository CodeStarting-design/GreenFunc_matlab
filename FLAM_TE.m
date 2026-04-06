function valid_poles = FLAM_TE(rho, h, er, freq)
    % --- 物理常数与波矢 ---
%     c = 3e8; 
    mu_0 = 4 * pi * 1e-7;
    eps_0 = 8.854e-12;
    omega = 2 * pi * freq;
    
    k0 = omega * sqrt(mu_0 * eps_0);
    k1 = omega * sqrt(mu_0 * er * eps_0);
    
%     if rho <= 0.05 * lambda
%         warning('距离过近：SDP-FLAM 适用于非近场区域 (rho > 0.05*lambda)。');
%     end


    % =====================================================
    % 步骤 1: FLAM - 寻找所有数学极点
    % =====================================================
%     disp('正在通过 FLAM 泰勒展开寻找复平面极点...');
    r2 = (k1^2 - k0^2) * h^2;
    N_max = 80; 
    a_points = -100:1:100; % 扫描实轴正半轴
    all_u_roots = [];
    
    for a = a_points
        % 构建局部泰勒展开多项式
        coeffs = zeros(1, N_max + 1);
        for n = 0:N_max
            coeffs(N_max + 1 - n) = cos(a + n * pi / 2) / factorial(n);
        end
        coeffs(N_max + 1) = coeffs(N_max + 1) + (a^2 / (2 * r2) - 1); 
        if N_max >= 1, coeffs(N_max) = coeffs(N_max) + a / r2; end   
        if N_max >= 2, coeffs(N_max - 1) = coeffs(N_max - 1) + 1 / (2 * r2); end 
        
        t_roots = roots(coeffs);
        u_roots = (t_roots + a) / 2;
        
        % 代回超越方程平方形式验证
        for i = 1:length(u_roots)
            u_t = u_roots(i);
            if abs(cos(2*u_t) + (2*u_t^2/r2 - 1)) < 1e-8
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
    % === 新增：k_rho 平面去重 ===
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
            
            % 2. 代回原始方程 Eq(5) 验算
            eq_top = 1i * kz0_top + kz1_val * cot(kz1_val * h);
            eq_bot = 1i * kz0_bot + kz1_val * cot(kz1_val * h);
            
            % 3. 分类判定
            if real(kp) > real(k0) && real(kp) < real(k1) && abs(imag(kp)) < 1e-4
                if abs(eq_top) < 1e-3  % 是真正的 SWP
                    SWP_list = [SWP_list; kp];
                end
            else
                if abs(eq_bot) < 1e-3  % 是真正的 LWP
                    LWP_list = [LWP_list; kp];
                end
            end
        end
    end
    
    fprintf('\n--- 物理极点识别结果 ---\n');
    fprintf('表面波极点 (SWP) 共 %d 个:\n', length(SWP_list));
    for i = 1:length(SWP_list)
        fprintf('  SWP %d: kp/k0 = %.4f + %.4fi\n', i, real(SWP_list(i))/real(k0), imag(SWP_list(i))/real(k0));
    end
    fprintf('漏波极点 (LWP) 共 %d 个\n', length(LWP_list));

    % 汇总物理极点
    physical_poles = [SWP_list; LWP_list];
    
    % --- 距离衰减判据 (洗刷远处的 LWP) ---
    valid_poles = physical_poles(abs(imag(physical_poles) .* rho) < 15);
%     fprintf('经过距离判据筛选后，最终参与留数计算的有效极点数: %d\n', length(valid_poles));
end