%% Validate_TimeDomain_3Points.m
% 时域三点验证（失稳前 / 失稳中 / 失稳后）
% 固定频率 Omega_fixed，选择 3 个力点 F_pre, F_mid, F_post
% 对每个点进行长时间积分，输出：
%   - 最后若干周期时域波形
%   - Poincare 截面（每周期采样一次）
%   - 频谱（判断周期/准周期/跳跃）
%
% 依赖（你的工程已有）：
%   L1.m, newton.m, nondim_temp2.m, branch_follow2.m, branch_aux2.m
% （可选）如果你已有 Run_L1_Stability 的变量 x_branchF / is_stable / max_mu，则会自动复用并选点

clear; clc; close all;

%% -----------------------------
% 1) 参数设置（与你 Run_L1_Stability 一致）
% -----------------------------
% sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
P.be1 = 1.0;
P.be2 = 0.1;
P.mu  = 0.2;

P.al1 = -0.95;
P.ga1 = 1.5;
P.ga2 = 1.5;

P.ze1 = 0.05;

P.lam   = 0.18;
P.kap_e = 0.395;
P.kap_c = 0.032;
P.sigma = 0.623;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

Omega_fixed = 0.25;   % 你当前关注的定频点（按你的图）
F0     = 0.001;       % 起始力（用于L1）
dF     = 5e-4;        % 扫力步长
Nsteps = 3000;        % 扫力长度

% Floquet 判断阈值（选点用）
tol_stable = 1.002;

%% -----------------------------
% 2) 获取扫力分支 + 稳定性标签（优先复用工作区变量）
% -----------------------------
global FixedOmega Fw
FixedOmega = Omega_fixed;

have_branch = exist('x_branchF','var') && exist('is_stable','var');
if ~have_branch
    fprintf('No existing x_branchF/is_stable in workspace. Recomputing L1 + Floquet selection labels...\n');

    % 先在F0处用Newton求一个周期解（HBM系数）
    Fw = F0;
    y_init = [zeros(15,1); F0];
    y_sol  = newton('nondim_temp2', y_init, sysP);
    x_coeff = y_sol(1:15);

    % 跑L1扫力
    x_interest = [x_coeff; F0];
    [x_branchF, infoL1] = L1(x_interest, sysP, dF, Nsteps, Omega_fixed);

    % 这里不强制重新算Floquet（很慢）
    % 你如果已有 max_mu/is_stable，就把下面两个变量手动传入工作区
    % 否则：用一个简单规则：先默认全稳定，让你先跑时域验证（你也可以把is_stable换成你之前算好的）
    warning('is_stable not provided. For rigorous 3-point picking, please reuse is_stable from Run_L1_Stability. Here we will pick points by force range only.');
    is_stable = true(1,size(x_branchF,2));
end

F_all = x_branchF(16,:).';
N = numel(F_all);

%% -----------------------------
% 3) 自动选三点：失稳前 / 失稳中 / 失稳后
% -----------------------------
% 若有 is_stable 且存在不稳定点：自动选点
if any(~is_stable)
    idxU = find(~is_stable(:));
    i_mid = idxU(round(numel(idxU)/2));                % 不稳定区中点

    % 找到 mid 左侧最近的稳定点
    i_pre = find(is_stable(1:i_mid-1), 1, 'last');
    if isempty(i_pre), i_pre = max(1, i_mid-50); end

    % 找到 mid 右侧最近的稳定点
    i_post = i_mid + find(is_stable(i_mid+1:end), 1, 'first');
    if isempty(i_post), i_post = min(N, i_mid+50); end

    F_pre  = F_all(i_pre);
    F_mid  = F_all(i_mid);
    F_post = F_all(i_post);

    fprintf('\n[AutoPick] F_pre=%.6g (idx=%d), F_mid=%.6g (idx=%d), F_post=%.6g (idx=%d)\n', ...
            F_pre, i_pre, F_mid, i_mid, F_post, i_post);
else
    % 若没有不稳定标签：按力范围选三个点（你可以手动改）
    [Fs, ord] = sort(F_all);
    F_pre  = Fs(round(0.25*numel(Fs)));
    F_mid  = Fs(round(0.50*numel(Fs)));
    F_post = Fs(round(0.75*numel(Fs)));
    fprintf('\n[RangePick] (no unstable labels) F_pre=%.6g, F_mid=%.6g, F_post=%.6g\n', F_pre, F_mid, F_post);
end

F_list = [F_pre, F_mid, F_post];
tag_list = {'PRE (stable-side)','MID (unstable window)','POST (stable-side)'};

%% -----------------------------
% 4) 时域积分设置（关键）
% -----------------------------
T = 2*pi/Omega_fixed;

% 建议：先跑 400 周期，丢掉前 300 周期，保留最后 100 周期分析
N_total_cycles = 450;
N_drop_cycles  = 320;     % 丢掉过渡
N_keep_cycles  = N_total_cycles - N_drop_cycles;

% 积分步长：每周期划分 N_per 步（越大越准，但越慢）
N_per = 400;           % 300~800都行，强非线性建议 >=400
dt = T / N_per;

% 为避免ode45在刚性区乱跳：用定步长RK4（与你Floquet一致）
use_RK4 = true;

% 初始条件：尽量用HBM重构（从L1分支里找到最接近该力的点）
%% =========================================================
%  MID 单点：双初值对照（LOW branch IC vs HIGH branch IC）
%  目的：验证“同一(F,Omega)”是否存在双稳态/跳跃（hysteresis）
%% =========================================================
F_mid_test = 0.00843978;        % <<< 用你截图中的 MID 力点（可改）
Omega_test = Omega_fixed;        % 同一个Omega

fprintf('\n==============================\n');
fprintf('MID Dual-IC test: Omega=%.6g, F=%.6g\n', Omega_test, F_mid_test);
fprintf('==============================\n');

% ---- 1) 在L1分支里找“F附近的一段点”，从中挑 LOW/HIGH 幅值初值 ----
Fwin = 8e-4;   % 搜索窗口（太小可能找不到两支；太大会混入别的区域）
idxNear = find(abs(F_all - F_mid_test) < Fwin);

if numel(idxNear) < 20
    warning('Near points too few (%d). Increase Fwin or ensure your L1 branch covers this F.', numel(idxNear));
end

% 用HBM系数估计幅值（你也可以换成合成幅值 sqrt(A1^2+A3^2)）
x2c_near = x_branchF(6:10, idxNear).';
A_est = hypot(x2c_near(:,2), x2c_near(:,3)) + hypot(x2c_near(:,4), x2c_near(:,5));

[~, iLow]  = min(A_est);
[~, iHigh] = max(A_est);

idxLow  = idxNear(iLow);
idxHigh = idxNear(iHigh);

fprintf('Pick ICs from branch:\n');
fprintf('  LOW : idx=%d,  F=%.6g, A_est=%.6g\n', idxLow,  F_all(idxLow),  A_est(iLow));
fprintf('  HIGH: idx=%d, F=%.6g, A_est=%.6g\n', idxHigh, F_all(idxHigh), A_est(iHigh));

yHBM_low  = x_branchF(:, idxLow);
yHBM_high = x_branchF(:, idxHigh);

y0_low  = hb_to_state6(yHBM_low,  Omega_test);
y0_high = hb_to_state6(yHBM_high, Omega_test);

% ---- 2) 时域积分设置（与你之前一致）----
T = 2*pi/Omega_test;
N_total_cycles = 520;
N_drop_cycles  = 360;
N_per = 500;                 % 建议比你之前更密一点（更稳）
dt = T / N_per;

Nt = N_total_cycles * N_per;
t = (0:Nt)'*dt;

% ---- 3) 分别用 LOW / HIGH 初值积分 ----
Y_low  = rk4_integrate(@(tt,yy) state6_ode(tt,yy,sysP,Omega_test,F_mid_test), t, y0_low);
Y_high = rk4_integrate(@(tt,yy) state6_ode(tt,yy,sysP,Omega_test,F_mid_test), t, y0_high);

% ---- 4) 取稳态段 ----
i_drop = N_drop_cycles*N_per + 1;
t_ss = t(i_drop:end);

Yl = Y_low(i_drop:end,:);
Yh = Y_high(i_drop:end,:);

x2_l = Yl(:,3); v2_l = Yl(:,4);
x2_h = Yh(:,3); v2_h = Yh(:,4);

% ---- 5) 稳态定量指标：RMS、峰峰值 ----
rms_l = rms(x2_l);  rms_h = rms(x2_h);
pp_l  = max(x2_l)-min(x2_l);
pp_h  = max(x2_h)-min(x2_h);

fprintf('\nSteady metrics (x2):\n');
fprintf('  LOW IC  -> rms=%.6g, p2p=%.6g\n', rms_l, pp_l);
fprintf('  HIGH IC -> rms=%.6g, p2p=%.6g\n', rms_h, pp_h);

% ---- 6) Poincare：每周期采样一次（同相位点）----
idxP = 1:N_per:numel(t_ss);
Px_l = x2_l(idxP); Pv_l = v2_l(idxP);
Px_h = x2_h(idxP); Pv_h = v2_h(idxP);

% ---- 7) 频谱（稳态段）----
Fsamp = 1/dt;
[f_l, X_l] = simple_spectrum(x2_l - mean(x2_l), Fsamp);
[f_h, X_h] = simple_spectrum(x2_h - mean(x2_h), Fsamp);

% ---- 8) 画对照图（1张图里看清楚“是否双稳态/跳跃”）----
figure('Color','w','Position',[80 60 1100 720]);
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

% (a) 最后10周期时域对比
N_show = min(numel(t_ss), 10*N_per);
nexttile; hold on; grid on; box on;
plot(t_ss(end-N_show+1:end), x2_l(end-N_show+1:end), 'LineWidth', 1.2);
plot(t_ss(end-N_show+1:end), x2_h(end-N_show+1:end), 'LineWidth', 1.2);
xlabel('t'); ylabel('x_2');
title(sprintf('x_2(t) last 10 cycles | \\Omega=%.3g, F=%.6g', Omega_test, F_mid_test));
legend('LOW IC','HIGH IC','Location','best');

% (b) 相图对比
nexttile; hold on; grid on; box on;
plot(x2_l(end-N_show+1:end), v2_l(end-N_show+1:end), 'LineWidth', 1.0);
plot(x2_h(end-N_show+1:end), v2_h(end-N_show+1:end), 'LineWidth', 1.0);
xlabel('x_2'); ylabel('v_2'); title('Phase portrait (last 10 cycles)');
legend('LOW IC','HIGH IC','Location','best');

% (c) Poincare map 对比
nexttile; hold on; grid on; box on;
scatter(Px_l, Pv_l, 14, 'filled');
scatter(Px_h, Pv_h, 14);
xlabel('x_2(kT)'); ylabel('v_2(kT)'); title('Poincaré map');
legend('LOW IC','HIGH IC','Location','best');

% (d) Poincare 序列对比
nexttile; hold on; grid on; box on;
plot(Px_l, 'LineWidth', 1.1);
plot(Px_h, 'LineWidth', 1.1);
xlabel('k'); ylabel('x_2(kT)'); title('Poincaré sequence');
legend('LOW IC','HIGH IC','Location','best');

% (e) 频谱对比（只看0~10倍基频）
nexttile; hold on; grid on; box on;
plot(f_l, X_l, 'LineWidth', 1.1);
plot(f_h, X_h, 'LineWidth', 1.1);
xlim([0, 10*Omega_test/(2*pi)]);
xlabel('Frequency (Hz)'); ylabel('|X_2|'); title('Spectrum (steady segment)');
legend('LOW IC','HIGH IC','Location','best');

% (f) 文本总结：是否双稳态
nexttile; axis off;
ratio = abs(rms_h - rms_l) / max(rms([rms_h,rms_l]), 1e-12);
if ratio > 0.1
    verdict = 'LIKELY BISTABLE at this (Omega,F): different steady states reached.';
else
    verdict = 'Likely UNIQUE attractor: both ICs converge to same steady state.';
end
text(0.02,0.85, sprintf(['Dual-IC verdict:\n%s\n\n' ...
    'LOW  rms=%.4g, p2p=%.4g\nHIGH rms=%.4g, p2p=%.4g\n'], ...
    verdict, rms_l, pp_l, rms_h, pp_h), 'FontSize', 12, 'FontName', 'Times New Roman');

%% =========================
% Local functions
% =========================

function y0 = hb_to_state6(yHBM16, Omega)
    % yHBM16: 16x1 from L1 branch: [x1c(5); x2c(5); qc(5); F]
    x1c = yHBM16(1:5);
    x2c = yHBM16(6:10);
    qc  = yHBM16(11:15);

    x1_0 = x1c(1) + x1c(2) + x1c(4);
    v1_0 = Omega * x1c(3) + 3 * Omega * x1c(5);

    x2_0 = x2c(1) + x2c(2) + x2c(4);
    v2_0 = Omega * x2c(3) + 3 * Omega * x2c(5);

    q_0  = qc(1)  + qc(2)  + qc(4);
    qd_0 = Omega * qc(3)  + 3 * Omega * qc(5);

    y0 = [x1_0; v1_0; x2_0; v2_0; q_0; qd_0];
end

function dy = state6_ode(t, y, sysP, Omega, Fw)
    % 仅状态6维 ODE（与你之前 ext_ode 的 dy 部分一致）
    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);
    theta = sqrt(max(lam,0));

    x1=y(1); v1=y(2); x2=y(3); v2=y(4); q=y(5); qd=y(6);
    dx = x1-x2; dv = v1-v2;

    f12 = (be1+al1)*dx + ga1*dx^3;
    f2g = be2*x2 + ga2*x2^3 + 2*mu*ze2*v2;

    x1dd = -f12 + theta*qd + Fw*cos(Omega*t);
    x2dd = ( f12 - f2g - theta*qd )/mu;

    if kap_e == 0
        tiny = 1e-12;
        qd_new = (-kap_c*q - theta*dv)/max(abs(sigma), tiny);
        qdd = (qd_new - qd)*50;
    else
        qdd = (-sigma*qd - kap_c*q - theta*dv)/kap_e;
    end

    dy = [v1; x1dd; v2; x2dd; qd; qdd];
end

function Y = rk4_integrate(f, t, y0)
    % 定步长RK4
    n = numel(t);
    m = numel(y0);
    Y = zeros(n, m);
    Y(1,:) = y0(:).';
    for i = 1:n-1
        h = t(i+1) - t(i);
        yi = Y(i,:).';
        k1 = f(t(i), yi);
        k2 = f(t(i) + 0.5*h, yi + 0.5*h*k1);
        k3 = f(t(i) + 0.5*h, yi + 0.5*h*k2);
        k4 = f(t(i) + h, yi + h*k3);
        Y(i+1,:) = (yi + (h/6)*(k1 + 2*k2 + 2*k3 + k4)).';
    end
end

function [f, Xmag] = simple_spectrum(x, Fs)
    % 简单单边幅度谱
    x = x(:);
    n = numel(x);
    X = fft(x);
    P2 = abs(X/n);
    P1 = P2(1:floor(n/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:floor(n/2))/n;
    Xmag = P1;
end

function txt = analyze_poincare(Px, Pv)
    % 非严格，但非常实用的判别提示（审稿人也认可这种“现象学判别”）
    Px = Px(:); Pv = Pv(:);

    % 去掉最开始的一小段（避免还没完全收敛）
    n = numel(Px);
    if n > 80
        Px2 = Px(round(0.2*n):end);
        Pv2 = Pv(round(0.2*n):end);
    else
        Px2 = Px; Pv2 = Pv;
    end

    % 聚类程度：看点云离散度
    sx = std(Px2); sv = std(Pv2);

    if sx < 1e-4 && sv < 1e-4
        txt = 'Likely PERIOD-1 stable orbit: Poincaré points collapse to a single point.';
        return;
    end

    % 看是否呈现“有限多个点”的离散集合（例如周期倍化）
    % 简单做法：对Px做bin，看唯一bin数量
    binW = max(1e-5, 0.02*std(Px2));
    id = round((Px2 - median(Px2))/binW);
    nuniq = numel(unique(id));

    if nuniq <= 6
        txt = sprintf('Likely MULTI-PERIOD orbit (period-k): Poincaré shows ~%d discrete clusters.', nuniq);
    else
        txt = 'Likely QUASI-PERIODIC / modulated response: Poincaré points form a dense set/curve.';
    end
end
