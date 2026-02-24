% =========================================================
% BG Model: Backward Sweep TF (add 1 electromech curve)
% 目标：复现 Wang(2017) 图(b) 的 4 条 L 曲线 + 额外加 1 条电路开启曲线
% =========================================================

clear; clc; close all;

% -------- Wang 图注给定参数（BG Model）--------
mu   = 0.2;     % m2/m1
beta = 2.0;     % 下层竖向线性刚度比 β
K1   = 1.0;     % 上层水平刚度比 K1
K2   = 0.2;     % 下层水平刚度比 K2
U    = 2.0;     % 几何尺度 U

% 曲线 1~4 的 L
L_list = [4/9, 0.5, 0.55, 0.6];

% -------- 用 LQZS=4/9 反推 v（使 alpha1(LQZS)=0）--------
% alpha1 = v - 2*K1*(1-L)/L
% 令 L=4/9 时 alpha1=0  =>  v = 2*K1*(1-L)/L = 2.5
v = 2.5;

% -------- 扫频设置（倒序）--------
global Fw FixedOmega
Fw = 0.005;         % Wang 图注：Fe = 0.005
FixedOmega = [];    % 扫频模式

Omega_Start = 10.0;
Omega_End   = 0.2;
Omega_Step  = -0.01;

nSteps = round((Omega_End - Omega_Start)/Omega_Step);  % 正数
nSteps = max(nSteps, 2000); % 保底

use_full13 = true;  % true: 1+3 合成；false: 仅基波

% =========================================================
% 额外第 5 条曲线：电路开启（你给的参数）
% 默认 L 用 4/9（你也可以改成 0.5/0.55/0.6）
% =========================================================
L_elec = 4/9;

cases = [];
for i = 1:numel(L_list)
    cases(end+1).name  = sprintf('Curve %d: L=%.4g (lam=0)', i, L_list(i));
    cases(end).L       = L_list(i);
    cases(end).lam     = 0.0;
    cases(end).kap_e   = 0.0;
    cases(end).kap_c   = 0.0;
    cases(end).sigma   = 0.0;
end

cases(end+1).name  = sprintf('Electromech: L=%.4g (lam=0.18,RLC)', L_elec);
cases(end).L       = L_elec;
cases(end).lam     = 0.18;
cases(end).kap_e   = 1.83;
cases(end).kap_c   = 0.01;
cases(end).sigma   = 0.54;

% =========================================================
% 开始循环计算并叠加绘图
% =========================================================
figure('Color','w'); ax = gca;
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'XScale','log');
xlim(ax, [Omega_End Omega_Start]);
yline(ax, 0, '-');
xlabel(ax, '\Omega (log scale)');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax, sprintf('BG Model: Backward Sweep TF (%d cases)', numel(cases)));

for ic = 1:numel(cases)
    L = cases(ic).L;

    % -------- Wang 等效线性/非线性系数 --------
    alpha1 = v    - 2*K1*(1-L)/L;
    alpha2 = beta - 2*K2*(1-L)/L;
    gamma1 = K1/(U^2 * L^3);
    gamma2 = K2/(U^2 * L^3);

    % -------- 映射到你的 sysP（关键对齐）--------
    P.be1 = 1.0;                 % 你的模型基准层间刚度
    P.al1 = alpha1 - P.be1;      % 保证 (be1 + al1) == alpha1
    P.be2 = alpha2;              % 对地线性项
    P.ga1 = gamma1;              % 上层三次
    P.ga2 = gamma2;              % 下层三次

    P.mu  = mu;
    P.ze1 = 0.05;                % Wang: ζ1=ζ2=0.05（你这里 ze1 对应 ζ2）

    % 电路参数（本 case）
    P.lam   = cases(ic).lam;
    P.kap_e = cases(ic).kap_e;
    P.kap_c = cases(ic).kap_c;
    P.sigma = cases(ic).sigma;

    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

    fprintf('\n[%d/%d] Running case: %s\n', ic, numel(cases), cases(ic).name);
    fprintf('   alpha1=%.6g, alpha2=%.6g, gamma1=%.6g, gamma2=%.6g\n', alpha1, alpha2, gamma1, gamma2);

    % -------- (1) 高频起点 --------
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;

    [x0_full, ok] = newton('nondim_temp2', y_init, sysP);
    if ~ok
        warning('Case "%s" 高频起点失败，跳过。', cases(ic).name);
        continue;
    end
    x0 = x0_full(1:15);

    % -------- (2) 第二点 --------
    Omega_Next = Omega_Start + Omega_Step;
    y_init2 = [x0; Omega_Next];

    [x1_full, ok] = newton('nondim_temp2', y_init2, sysP);
    if ~ok
        warning('Case "%s" 第二点失败，跳过。', cases(ic).name);
        continue;
    end
    x1 = x1_full(1:15);

    % -------- (3) 弧长延拓 --------
    [x_res, ~] = branch_follow2('nondim_temp2', nSteps, Omega_Start, Omega_Next, x0, x1, sysP);

    % -------- (4) 严格按 ft 定义算 TF --------
    Om  = x_res(16,:).';           % Nx1
    be2 = sysP(2);
    mu2 = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);

    x2 = x_res(6:10,:).';          % Nx5

    % x2' HB系数
    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);

    % x2^3 HB系数 (AFT)
    x2_cub = cubic_proj_013_batch(x2);

    % ft HB系数
    ft = be2*x2 + ga2*x2_cub + 2*mu2*ze2*x2_dot;  % Nx5

    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));

    if use_full13
        ft_amp = hypot(ft1, ft3);
    else
        ft_amp = ft1;
    end

    TF    = ft_amp ./ Fw;
    TF_dB = 20*log10(max(TF, 1e-300));

    okp = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om = Om(okp); TF_dB = TF_dB(okp);

    % ✅用线连接（你现在已经不需要 scatter 了）
    plot(ax, Om, TF_dB, '-', 'LineWidth', 1.8, 'DisplayName', cases(ic).name);
end

legend(ax,'Location','northeast');

% ============ AFT 批处理（和 nondim_temp2 一致）===========
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
