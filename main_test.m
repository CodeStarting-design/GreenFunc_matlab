% =========================================================================
% 高阶索末菲积分算法: 矩形绕道 + 双指数 Tanh-Sinh 分区外推法
% 精确计算无损耗分层介质空间域格林函数 GAxx 与 Gphi
% =========================================================================
% clear; clc;

% --- 1. 物理参数设置 ---
f = 8e9;                  % 8 GHz
c = 3e8;
lambda0 = c / f;
h = lambda0;
% rho = 5*lambda0;      % 观测距离
rhomin=0.01*lambda0;
rhomax=10*lambda0;
Nx=1000;
rho=logspace(log10(rhomin), log10(rhomax), Nx);

% 无需人工损耗，采用纯实数介电常数
er = 4;           

% fprintf('==================================================\n');
% fprintf('  高阶索末菲积分计算 (矩形绕道 + Levin-Sidi 外推)\n');
% fprintf('  频率: %g GHz, 距离 rho: %g m\n', f/1e9, rho);
% fprintf('  介电常数 er = %g (无损耗)\n', er);
% fprintf('==================================================\n');

% --- 2. 波矢计算 ---
mu_0 = 4 * pi * 1e-7; 
eps_0 = 8.854e-12; 
omega = 2 * pi * f;
k0 = omega * sqrt(mu_0 * eps_0);
k1 = omega * sqrt(mu_0 * er * eps_0);

% 设定实轴绕道的切分点 a
% 必须保证 a 大于所有的物理极点 (SWP 通常在 k0 到 k1 之间)
a_point = 1.2 * real(k1); 

% =========================================================================
% 3. 计算 GAxx
% =========================================================================
% disp('>>> 正在计算 GAxx ...');
valid_poles1 = FLAM_TE(rho(1), h, er, f);
valid_poles2 = FLAM_TM(rho(1), h, er, f);
valid_poles=[valid_poles1;valid_poles2];
G_A = zeros(Nx,1);
G_A_DCIM = zeros(Nx,1);
for i=1:length(rho)

% integrand_GA = @(kr) Integrand_GAxx(kr, rho(i), k0, k1, h);
% 
% % [0, a] 区间: 复平面矩形绕道 (避开实轴极点)
% GA_near = IntegrateSpectralNearField(integrand_GA, a_point, k1);
% 
% % [a, inf) 区间: 贝塞尔零点分区 + Tanh-Sinh + Levin-Sidi 外推
% GA_far  = IntegrateSpectralFarField(integrand_GA, a_point, rho(i));
% 
% GA_total(i) = GA_near + GA_far;

G_A(i) = calculate_GAxx_SDP_FLAM(valid_poles1, rho(i), h, er, f);

G_A_DCIM(i) = calculate_GAxx_DCIM2(valid_poles1, rho(i), h, er, f);
% % =========================================================================
% % 4. 计算 Gphi
% % =========================================================================
% % disp('>>> 正在计算 Gphi ...');
% integrand_Gphi = @(kr) Integrand_Gphi(kr, rho(i), k0, k1, h, er);
% 
% % [0, a] 区间: 复平面矩形绕道
% Gphi_near = IntegrateSpectralNearField(integrand_Gphi, a_point, k1);
% 
% % [a, inf) 区间: 贝塞尔零点分区 + Levin-Sidi 外推
% Gphi_far  = IntegrateSpectralFarField(integrand_Gphi, a_point, rho(i));
% 
% Gphi_total(i) = Gphi_near + Gphi_far;
% % fprintf('  Gphi = %e + %ei\n', real(Gphi_total), imag(Gphi_total));
% % fprintf('  |Gphi| = %e\n', abs(Gphi_total));
% % fprintf('==================================================\n');
% 
% 
% G_phi(i) = calculate_Gphi_SDP_FLAM(valid_poles, rho(i), h, er, f);
end

figure(1)
% p1=loglog(rho/lambda0,abs(GA_total),'k-','LineWidth', 2);
% hold on
p2=loglog(rho/lambda0,abs(G_A)/mu_0,'r--','LineWidth', 2);
hold on
p3=loglog(rho/lambda0,abs(G_A_DCIM)/mu_0,'m--','LineWidth', 2);

% figure(2)
% p4=loglog(rho/lambda0,abs(Gphi_total)/eps_0,'k-','LineWidth', 2);
% hold on
% p5=loglog(rho/lambda0,abs(G_phi),'r--','LineWidth', 2);

%% =========================================================================
% 积分核函数定义
% =========================================================================

% --- GAxx 的积分核 ---
function val = Integrand_GAxx(kr, rho, k0, k1, h)
    kz0 = sqrt(k0^2 - kr.^2);
    kz0(imag(kz0) > 0) = -kz0(imag(kz0) > 0); % 强制衰减条件
    kz1 = sqrt(k1^2 - kr.^2);
    
    f1 = kz1 .* cos(kz1 .* h) + 1i .* kz0 .* sin(kz1 .* h);
    K11 = sin(kz1 .* h) ./ f1;
    val = (1 / (2*pi)) .* K11 .* besselj(0, kr .* rho) .* kr;
end

% --- Gphi 的积分核 ---
function val = Integrand_Gphi(kr, rho, k0, k1, h, er)
    kz0 = sqrt(k0^2 - kr.^2);
    kz0(imag(kz0) > 0) = -kz0(imag(kz0) > 0); 
    kz1 = sqrt(k1^2 - kr.^2);
    
    f1 = kz1 .* cos(kz1 .* h) + 1i .* kz0 .* sin(kz1 .* h);
    f2 = er .* 1i .* kz0 .* cos(kz1 .* h) - kz1 .* sin(kz1 .* h);
    
    num = 1i .* kz0 .* sin(kz1 .* h) .* cos(kz1 .* h) - kz1 .* (sin(kz1 .* h)).^2;
    K22 = num ./ (f1 .* f2);
    val = (1 / (2*pi)) .* K22 .* besselj(0, kr .* rho) .* kr;
end

%% =========================================================================
% 积分调度算法 (近场绕道 + 远场分区外推)
% =========================================================================

function val = IntegrateSpectralNearField(f, a, k1)
    % 高度 h_detour 设置为介质波数的 1/1000
    h_detour = 1.0e-3 * real(k1);
    tol = 1e-6;
    
    % Path 1: (0 -> i*h)
    p1 = integral(f, 0, 1i*h_detour, 'AbsTol', tol, 'RelTol', tol, 'ArrayValued', true);
    % Path 2: (i*h -> a + i*h)
    p2 = integral(f, 1i*h_detour, a + 1i*h_detour, 'AbsTol', tol, 'RelTol', tol, 'ArrayValued', true);
    % Path 3: (a + i*h -> a)
    p3 = integral(f, a + 1i*h_detour, a, 'AbsTol', tol, 'RelTol', tol, 'ArrayValued', true);
    
    val = p1 + p2 + p3;
end

function val = IntegrateSpectralFarField(f, a, rho)
    tol = 1e-6;
    
    if rho < 1e-15
        val = integral(f, a, Inf, 'RelTol', tol, 'AbsTol', 1e-12, 'ArrayValued', true);
        return;
    end
    
    % 动态寻找下一个贝塞尔零点
    b = BesselJ_NextZero(a, rho);
    q = pi / rho; % 半周期长度
    
    if b > a
        bridge = integral(f, a, b, 'RelTol', tol, 'AbsTol', 1e-12, 'ArrayValued', true);
    else
        bridge = 0;
    end
    
    % 执行分区外推
    tail = PartExtrap(f, b, q, tol);
    val = bridge + tail;
end

% --- 动态贝塞尔零点搜索 (代替硬编码表) ---
function b = BesselJ_NextZero(a, rho)
    z_a = a * rho;
    m = ceil(z_a / pi + 0.25);
    z_guess = (m - 0.25) * pi;
    if z_guess <= z_a
        m = m + 1;
        z_guess = (m - 0.25) * pi;
    end
    z_root = fzero(@(z) besselj(0, z), z_guess);
    b = z_root / rho;
end

%% =========================================================================
% Tanh-Sinh 双指数积分与 Levin-Sidi 外推核心算法
% =========================================================================

function val = PartExtrap(f, a, q, tol)
    kmax = 10;
    X = zeros(kmax + 2, 1); X(1) = a;
    A = complex(zeros(kmax + 1, 1)); B = complex(zeros(kmax + 1, 1));
    s = 0; old_val = 1e300; val = 0;
    
    for kk = 1:(kmax + 1)
        X(kk+1) = X(kk) + q;
        u = TanhSinh(f, X(kk), X(kk+1), tol);
        if abs(u) == 0, continue; end
        
        s = s + u;
        w = u;
        [val, A, B] = LevinSidi(kk, s, w, X, A, B);
        
        if kk > 2 && abs(val - old_val) < tol * abs(val)
            break;
        end
        old_val = val;
    end
end

function [result, A, B] = LevinSidi(k, Sk, wk, X, A, B)
    B(k) = 1.0 / wk;
    A(k) = Sk * B(k);
    for j = 1:(k-1)
        d = 1.0 / X(k+1) - 1.0 / X(k-j+1);
        A(k-j) = (A(k-j+1) - A(k-j)) / d;
        B(k-j) = (B(k-j+1) - B(k-j)) / d;
    end
    result = A(1) / B(1);
end

function val = TanhSinh(f0, a, b, tol)
    f = @(c, d) f0(c + d);
    eta = 1.0; maxlev = 5;
    sigma = (b - a) / 2.0; gamma = (b + a) / 2.0;
    
    s = eta * f(gamma, 0.0);
    h = 0.5; eh = exp(h);
    
    [n, s] = TruncIndex(f, eh, s, sigma, eta, a, b);
    old_val = sigma * h * s; val = old_val;

    for m = 1:maxlev
        e2h = eh; h = h / 2.0; eh = exp(h);
        s_new = PartSum(f, eh, e2h, n, sigma, eta, a, b);
        val = old_val / 2.0 + sigma * h * s_new;
        if abs(val - old_val) <= tol * abs(val), break; end
        old_val = val; n = n * 2;
    end
end

function [n, s] = TruncIndex(f, eh, s, sigma, eta, a, b)
    nmax = 24; kappa = 1.0e-15; ekh = eh; n = 0;
    for i = 0:nmax
        n = i;
        t = Term(f, ekh, sigma, eta, a, b);
        s = s + t;
        if abs(t) <= kappa * abs(s), break; end
        ekh = ekh * eh;
    end
    n = n - 1;
end

function s = PartSum(f, eh, e2h, n, sigma, eta, a, b)
    ekh = eh;
    s = Term(f, ekh, sigma, eta, a, b);
    for k = 2:n
        ekh = ekh * e2h;
        t = Term(f, ekh, sigma, eta, a, b);
        s = s + t;
    end
end

function t = Term(f, ekh, sigma, eta, a, b)
    inv_ekh = 1.0 / ekh;
    q = exp(-eta * (ekh - inv_ekh));
    delta = 2.0 * q / (1.0 + q);
    w = eta * (ekh + inv_ekh) * delta / (1.0 + q);
    t = w * (f(a, sigma * delta) + f(b, -sigma * delta));
end