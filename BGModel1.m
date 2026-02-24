% =========================================================
% 目标：Wang(2017) BG 模型：改变 L 画 4 条 TF 曲线
%      (σ=0.2, U=2, β=2, K1=1, K2=0.2, ζ1=ζ2=0.05, Fe=0.005)
%
% 说明：
%   - 按图注：Curve1 L=4/9; Curve2 L=0.5; Curve3 L=0.55; Curve4 L=0.6
%   - v 取 2.5（由 LQZS=4/9 令 alpha1=0 反推得到），并对四条曲线保持不变
%   - 你的 sysP 对齐规则：
%       (be1+al1) = alpha1(Wang)  => al1 = alpha1 - be1
%       be2       = alpha2(Wang)
%       ga1,ga2   = gamma1,gamma2(Wang)
% =========================================================

clear; clc; close all;

% -------- Wang 图注给定参数（BG Model）--------
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比（Wang 的 β）
K1   = 1.0;     % 上层水平弹簧刚度比（Wang 的 K1）
K2   = 0.2;     % 下层水平弹簧刚度比（Wang 的 K2）
U    = 2.0;     % 几何非线性尺度参数
zeta = 0.05;    % 图注 ζ1=ζ2=0.05
Fw0  = 0.005;   % 图注 Fe=0.005（你代码里用全局 Fw）

% -------- 关键：v 固定（由 LQZS=4/9 使 alpha1=0 反推）--------
% alpha1 = v - 2*K1*(1-L)/L
% 令 L=4/9 时 alpha1=0 => v = 2*K1*(1-L)/L = 2.5
v = 2.5;

% -------- 四个 L --------
L_list = [4/9, 0.5, 0.55, 0.6];

% -------- 扫频设置：倒序扫频 --------
Omega_Start = 10.0;
Omega_End   = 0.2;
Omega_Step  = -0.01;
nStepsArc   = 5000;

% -------- 全局（你的 FRF/HBM 框架）--------
global Fw FixedOmega
Fw = Fw0;
FixedOmega = [];   % 扫频模式

% ======== 画图窗口（四条曲线叠加）========
figure('Color','w');
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'XScale','log');
xlim(ax, [Omega_End Omega_Start]);
yline(ax, 0, '-');
xlabel(ax, '\Omega (log scale)');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax, 'BG Model: Backward Sweep TF (4 L cases)');

use_full13 = true;  % true=1+3 合成；false=仅基波

% ======== 主循环：对每个 L 跑一次 continuation ========
for i = 1:numel(L_list)
    L = L_list(i);

    % ---- Wang 的等效线性/非线性系数（随 L 变化）----
    alpha1 = v    - 2*K1*(1-L)/L;     % 上层层间线性系数（乘 (x1-x2)）
    alpha2 = beta - 2*K2*(1-L)/L;     % 下层对地线性系数（乘 x2）
    gamma1 = K1/(U^2 * L^3);          % 上层三次
    gamma2 = K2/(U^2 * L^3);          % 下层三次

    % ---- 塞到你的 sysP（严格保持你的代码结构）----
    P.be1 = 1.0;
    P.al1 = alpha1 - P.be1;           % ★确保 (be1+al1)=alpha1
    P.be2 = alpha2;
    P.ga1 = gamma1;
    P.ga2 = gamma2;

    P.mu  = mu;
    P.ze1 = zeta;

    % 电路关闭
    P.lam   = 0.0;
    P.kap_e = 0.0;
    P.kap_c = 0.0;
    P.sigma = 0.0;

    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

    fprintf('\n===== Case %d: L = %.6f =====\n', i, L);
    fprintf('alpha1=%.6f, alpha2=%.6f, gamma1=%.6f, gamma2=%.6f\n', alpha1, alpha2, gamma1, gamma2);

    % ---- 1) 高频起点（Newton）----
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;  % [15 HB + Omega]

    [x0_full, ok] = newton('nondim_temp2', y_init, sysP);
    if ~ok, error('Case %d: 高频起点求解失败！', i); end
    x0 = x0_full(1:15);

    % ---- 2) 第二点 ----
    Omega_Next = Omega_Start + Omega_Step;
    y_init2 = [x0; Omega_Next];

    [x1_full, ok] = newton('nondim_temp2', y_init2, sysP);
    if ~ok, error('Case %d: 第二个点求解失败！', i); end
    x1 = x1_full(1:15);

    % ---- 3) 弧长延拓（倒序扫频）----
    [x_res, ~] = branch_follow2('nondim_temp2', nStepsArc, Omega_Start, Omega_Next, x0, x1, sysP);

    % ---- 4) 计算 TF（严格 ft=be2*x2 + ga2*x2^3 + 2µζ2*x2'）----
    Om  = x_res(16,:).';      % Nx1
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    x2 = x_res(6:10,:).';     % Nx5

    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);

    x2_cub = cubic_proj_013_batch(x2);

    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;  % Nx5

    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));

    if use_full13
        ft_amp = hypot(ft1, ft3);
    else
        ft_amp = ft1;
    end

    TF_dB = 20*log10(max(ft_amp./Fw, 1e-300));

    ok = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om = Om(ok);
    TF_dB = TF_dB(ok);

    % ---- 5) 画线（你说要“直线连接”）----
    % 注意：这里不 sort，保持延拓输出顺序（倒序扫频时一般是连续的）
    plot(ax, Om, TF_dB, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('Curve %d: L=%.4g', i, L));
end

legend(ax, 'Location', 'best');


%% ============ AFT：批量计算 x^3 的 HB 系数（0/1/3） ============
function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    X_time  = (T_mat * U.').';       % N x Nt
    X3_time = X_time.^3;
    cubic   = (T_inv * X3_time.').'; % N x 5
end

function [N, T_mat, T_inv] = get_AFT_matrices_local()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64;
        t = (0:pN-1)'*(2*pi/pN);
        c1=cos(t); s1=sin(t); c3=cos(3*t); s3=sin(3*t); dc=ones(pN,1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end
    N = pN; T_mat = pT; T_inv = pTinv;
end
