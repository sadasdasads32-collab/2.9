%% Verify_ODE45_vs_FRF.m
% ODE45（时域） vs FRF-HBM（频域）同频点验证（适配你当前代码）
%
% sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
% FRF(sysP) 输出 16xN:
%   1:5   x1=[x10 a11 b11 a13 b13]
%   6:10  x2=[x20 a21 b21 a23 b23]
%   11:15 q =[q0  aq1 bq1 aq3 bq3]
%   16    Omega
%
% ODE 状态：y = [x1; v1; x2; v2; q; qd]
%
% 依赖：FRF.m, branch_follow2.m, branch_aux2.m, newton.m, nondim_temp2.m

clc; close all;

%% ===================== 0) 工况设置 =====================
global Fw
Fw = 0.005;

% ---- 这里填你“当前能跑通”的参数 ----
P.be1   = 1.0;
P.be2   = 0.1;
P.mu    = 0.2;
P.al1   = -0.95;

% QZS（示例，按你自己的来）
k1=1; k2=0.8; L=4/9; U=2;
P.ga1   = k1/(U^2*L^3);
P.ga2   = k2/(U^2*L^3);

P.ze1   = 0.05;   % 代码里叫 ze1，但物理上是 ζ2（下层对地阻尼）

% 电路（示例：你可以换成你现在要验证的那组）
P.lam   = 0.18;
P.kap_e = 0.2;
P.kap_c = 0.5;
P.sigma = 1.0;

sysP = [P.be1; P.be2; P.mu; P.al1; P.ga1; P.ze1; ...
        P.lam; P.kap_e; P.kap_c; P.sigma; P.ga2];

% -------- 选择验证频率点（建议选峰附近）--------
Omega_test = 0.80;

% 耦合符号（若发现ODE/FRF明显反相或误差巨大，就改 -1 再试）
sgn_em = +1;

%% ===================== 1) ODE45 配置 =====================
T_periods = 400;          % 总周期数（强非线性可 300~800）
N_last    = 30;           % 取最后 N_last 周期作为稳态
T1 = 2*pi/Omega_test;
tspan = [0, T_periods*T1];

% 初值（从静止启动；若收敛慢给个小扰动）
y0 = zeros(6,1);
% y0(1) = 1e-4;

opts = odeset('RelTol',1e-7,'AbsTol',1e-9, 'MaxStep', T1/50);

fprintf('====================================================\n');
fprintf('ODE45 vs FRF verification\n');
fprintf('Omega_test=%.4f, Fw=%.4f, sgn_em=%+d\n', Omega_test, Fw, sgn_em);
fprintf('lam=%.4f kap_e=%.4f kap_c=%.4f sigma=%.4f\n', P.lam, P.kap_e, P.kap_c, P.sigma);
fprintf('====================================================\n');

%% ===================== 2) 跑 ODE45 =====================
fprintf('Running ODE45...\n');
[t,y] = ode45(@(t,y) odesys_nd(t,y,sysP,Omega_test,Fw,sgn_em), tspan, y0, opts);

% 截取稳态
I = find(t > t(end) - N_last*T1);
t_s  = t(I);
x1_s = y(I,1);
v1_s = y(I,2);

%% ===================== 3) ODE 幅值提取（A_pp 与 A_1x） =====================
A_pp = (max(x1_s)-min(x1_s))/2;

x = x1_s - mean(x1_s);
C = 2/numel(t_s) * sum(x .* cos(Omega_test*t_s));
S = 2/numel(t_s) * sum(x .* sin(Omega_test*t_s));
A_1x = hypot(C,S);

fprintf('\n[ODE steady]\n');
fprintf('A_pp = %.6g\n', A_pp);
fprintf('A_1x = %.6g (sync demod, 1x)\n', A_1x);

%% ===================== 4) 跑 FRF 并在 Omega_test 取幅值 =====================
fprintf('\nRunning FRF(sysP)...\n');
x_frf = FRF(sysP);

FRF_ok = ~isempty(x_frf) && size(x_frf,1)>=16 && size(x_frf,2)>=3;
A_frf_1x = NaN; err_pct = NaN;

if FRF_ok
    Om = x_frf(16,:).';
    Amp1x = sqrt(x_frf(2,:).^2 + x_frf(3,:).^2).'; % x1 1x 幅值

    % 排序（防倒序）
    [Om, idx] = sort(Om);
    Amp1x = Amp1x(idx);

    % 检查 FRF 是否覆盖 Omega_test
    if Omega_test < min(Om) || Omega_test > max(Om)
        fprintf('[FRF] WARNING: Omega_test=%.4f 不在 FRF 范围 [%.4f, %.4f] 内，不做插值。\n', ...
            Omega_test, min(Om), max(Om));
    else
        A_frf_1x = interp1(Om, Amp1x, Omega_test, 'linear');
        err_pct = abs(A_1x - A_frf_1x)/max(abs(A_frf_1x),1e-12)*100;

        fprintf('\n[FRF @ Omega_test]\n');
        fprintf('A_FRF_1x = %.6g\n', A_frf_1x);
        fprintf('\n[Error]\n');
        fprintf('Relative error = %.2f %%\n', err_pct);
    end
else
    fprintf('[FRF] FAILED or too few points. 请先确保 FRF 能扫过 Omega_test。\n');
end

%% ===================== 5) 图：稳态时域 / 相图 =====================
figure('Color','w','Name','Steady time response');
plot(t_s, x1_s, 'LineWidth',1.2); grid on;
xlabel('t'); ylabel('x_1');
title(sprintf('ODE steady x_1(t), \\Omega=%.4f', Omega_test));

figure('Color','w','Name','Phase portrait');
plot(x1_s, v1_s, 'LineWidth',1.0); grid on;
xlabel('x_1'); ylabel('x_1''');
title('Phase portrait (check periodicity)');

%% ===================== 6) 频谱（看谐波含量） =====================
figure('Color','w','Name','Spectrum (x1 steady)');
xw = x - mean(x);
N = numel(xw);
dt = mean(diff(t_s));
fs = 1/dt;

X = fft(xw .* hann(N));
f = (0:N-1)*(fs/N);

f_max = min(max(f), 10*(Omega_test/(2*pi)));
Iplot = (f <= f_max);

plot(f(Iplot), 20*log10(abs(X(Iplot))/max(abs(X(Iplot))+eps) + eps), 'LineWidth',1.2);
grid on;
xlabel('Frequency (Hz)'); ylabel('Normalized magnitude (dB)');
title('Spectrum (normalized)');

%% ===================== 7) FRF 曲线 + 验证点标记（可选） =====================
if FRF_ok
    figure('Color','w','Name','FRF and verification point');
    dB = 20*log10(Amp1x + eps);
    plot(Om, dB, 'LineWidth',1.4); grid on; hold on;
    xline(Omega_test,'r--','LineWidth',1.3);

    if isfinite(A_frf_1x) && isreal(A_frf_1x) && A_frf_1x>0
        yline(20*log10(A_frf_1x+eps),'r:','LineWidth',1.3);
    end
    xlabel('\Omega'); ylabel('20log_{10}(|X_1|) (dB)');
    title('FRF (HBM) with ODE verification point');
end

%% =========================================================
% 时域方程：与你 nondim_temp2 结构对齐
% =========================================================
function dydt = odesys_nd(t, y, sysP, Omega, Fw, sgn_em)
    % y = [x1; v1; x2; v2; q; qd]
    x1=y(1); v1=y(2);
    x2=y(3); v2=y(4);
    q =y(5); qd=y(6);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5);
    ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10);
    ga2=sysP(11);

    dx  = x1 - x2;
    f12 = (be1 + al1)*dx + ga1*dx^3;
    f2g = be2*x2 + ga2*x2^3 + 2*mu*ze2*v2;

    th = sqrt(max(lam,0));

    % 机械
    x1dd = -f12 + sgn_em*(th*qd) + Fw*cos(Omega*t);
    x2dd = ( f12 - f2g - sgn_em*(th*qd) ) / mu;

    % 电路：kap_e*qdd + sigma*qd + kap_c*q = -th*(v1 - v2)
    if kap_e == 0
        tiny = 1e-12;
        qd_new = (-kap_c*q - th*(v1 - v2)) / max(abs(sigma), tiny);
        qdd = (qd_new - qd) * 50;
    else
        qdd = (-sigma*qd - kap_c*q - th*(v1 - v2)) / kap_e;
    end

    dydt = [v1; x1dd; v2; x2dd; qd; qdd];
end
