%% Verify_Point.m - 审稿人级验证（ODE45 vs HBM-FRF 同频点对比）
% 功能：
% 1) 固定一个工况(sysP)和一个频率 Omega_test
% 2) 用 ODE45 跑长时程到稳态，截取最后 N 个周期
% 3) 提取稳态幅值：
%    - A_pp = (max-min)/2
%    - A_1x = 1×Ω 同步检波幅值（推荐，与 HBM 1次谐波幅值同口径）
% 4) 再调用 FRF(sysP)，用插值取 Omega_test 处的 FRF 1×幅值
% 5) 输出误差百分比，画时域、相图、频谱
%
% 依赖：
% - odesys_phys(t,y,sysP,Omega,Fw)  (你的时域方程，y至少包含 x1 与 x1')
% - FRF(sysP) -> 16xN，其中第16行 Omega；第2/3行为 x1 1次谐波 cos/sin 系数
%
% 注意：
% - 若你的状态顺序里 y(:,2) 不是 x1'，请改 IDX_V1
%
% Author: ChatGPT
% Date: 2026-02-10

clc; close all;

%% -----------------------------
% 0) 与扫频一致的工况参数
% -----------------------------
global Fw
Fw = 0.05;

% sysP 顺序：
% [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
P.be1   = 1.0;
P.be2   = 1.0;
P.mu    = 0.2;
P.al1   = 0.1;
P.ga1   = 0.1;
P.ze1   = 0.05;
P.lam   = 0.5;      % 选择你出现峰分裂的 λ
P.kap_e = 2.04;
P.kap_c = 1.0;
P.sigma = 0.5;
P.ga2   = 0.2;

sysP_verify = [P.be1; P.be2; P.mu; P.al1; P.ga1; P.ze1; P.lam; ...
               P.kap_e; P.kap_c; P.sigma; P.ga2];

%% -----------------------------
% 1) 选择要验证的频率点（建议选峰值频率）
% -----------------------------
% 你之前输出过 λ=0.5 的分裂峰：0.5817 和 0.7716
Omega_test = 0.7716;   % <-- 改这里即可，例如 0.5817

%% -----------------------------
% 2) ODE45 时域仿真设置
% -----------------------------
T_periods = 350;                 % 跑的周期数（强耦合/分裂建议 300~600）
T_end = T_periods * (2*pi/Omega_test);
tspan = [0, T_end];

% 初值：从静止启动（如收敛慢可给小扰动）
y0 = zeros(6,1);

options = odeset('RelTol',1e-7, 'AbsTol',1e-9);

fprintf('====================================================\n');
fprintf('ODE45 verification: lambda=%.2f, Omega=%.4f, Fw=%.4f\n', P.lam, Omega_test, Fw);
fprintf('Params: kap_e=%.3f kap_c=%.3f sigma=%.3f\n', P.kap_e, P.kap_c, P.sigma);
fprintf('====================================================\n');

fprintf('Running ODE45... (this may take a while)\n');
[t, y] = ode45(@(t,y) odesys_phys(t,y,sysP_verify,Omega_test,Fw), tspan, y0, options);

%% -----------------------------
% 3) 截取最后 N 个周期作为稳态
% -----------------------------
N_last = 30;   % 最后 20~50 周期
T1 = 2*pi/Omega_test;
I_steady = find(t > t(end) - N_last*T1);

t_s  = t(I_steady);
x1_s = y(I_steady, 1);

% 速度列号（默认 y(:,2)=x1'；若不对，请改）
IDX_V1 = 2;
v1_s = y(I_steady, IDX_V1);

%% -----------------------------
% 4) 幅值计算（pp/2 与 同步检波 1×）
% -----------------------------
% (a) 峰峰/2 幅值
A_pp = (max(x1_s) - min(x1_s))/2;
dB_pp = 20*log10(A_pp + eps);

% (b) 同步检波提取 1×Ω 幅值（与 HBM 同口径）
x = x1_s - mean(x1_s);
w = Omega_test;

C = 2/numel(t_s) * sum( x .* cos(w*t_s) );
S = 2/numel(t_s) * sum( x .* sin(w*t_s) );
A_1x = sqrt(C^2 + S^2);
dB_1x = 20*log10(A_1x + eps);

%% -----------------------------
% 5) 运行 FRF 并插值取同频点幅值（1×谐波）
% -----------------------------
fprintf('\nRunning FRF(...) to get frequency-domain amplitude...\n');

x_frf = [];
try
    x_frf = FRF(sysP_verify);         % 先试列向量
catch
    try
        x_frf = FRF(sysP_verify.');   % 再试行向量
    catch ME
        warning('FRF failed: %s', ME.message);
    end
end

FRF_ok = ~isempty(x_frf) && size(x_frf,1) >= 16;

if FRF_ok
    Omega_frf = x_frf(16,:);
    AmpX1_frf = sqrt(x_frf(2,:).^2 + x_frf(3,:).^2);          % 与你扫频一致
    dB_frf    = 20*log10(AmpX1_frf + eps);

    A_frf_1x  = interp1(Omega_frf, AmpX1_frf, Omega_test, 'linear', 'extrap');
    dB_frf_1x = 20*log10(A_frf_1x + eps);

    err_pct = abs(A_1x - A_frf_1x) / max(A_frf_1x, 1e-12) * 100;
else
    A_frf_1x = NaN; dB_frf_1x = NaN; err_pct = NaN;
end

%% -----------------------------
% 6) 打印结果（核心结论）
% -----------------------------
fprintf('\n================== Verification Result ==================\n');
fprintf('lambda = %.2f, Omega_test = %.4f\n', P.lam, Omega_test);

fprintf('\n[ODE45 steady-state]\n');
fprintf('A_pp  = %.6g  => dB = %.4f dB\n', A_pp, dB_pp);
fprintf('A_1x  = %.6g  => dB = %.4f dB   (recommended, 1× amplitude)\n', A_1x, dB_1x);

if FRF_ok
    fprintf('\n[FRF (HBM) @ Omega_test (interpolated)]\n');
    fprintf('A_FRF = %.6g  => dB = %.4f dB\n', A_frf_1x, dB_frf_1x);

    fprintf('\n[Error]\n');
    fprintf('Relative error (ODE 1× vs FRF) = %.2f %%\n', err_pct);
    if err_pct < 5
        fprintf('PASS: error < 5%%, frequency-domain result is credible.\n');
    else
        fprintf('NOTE: error >= 5%%, check harmonics/multi-solution/ODE convergence.\n');
    end
else
    fprintf('\n[FRF failed] => cannot compute numeric error automatically.\n');
end
fprintf('=========================================================\n');

%% -----------------------------
% 7) 绘图：稳态时域 + 幅值线
% -----------------------------
figure('Name','Steady-State Time Response','Color','w');
plot(t_s, x1_s, 'b', 'LineWidth', 1.2); grid on; hold on;
yline(A_pp,  'r--', 'LineWidth', 1.5);
yline(-A_pp, 'r--', 'LineWidth', 1.5);
title(sprintf('Steady-State x_1(t) @ \\lambda=%.2f, \\Omega=%.4f', P.lam, Omega_test));
xlabel('t'); ylabel('x_1');
legend('x_1(t)','\pm A_{pp}','Location','best');

%% -----------------------------
% 8) 相图：检查周期性
% -----------------------------
figure('Name','Phase Portrait (Check periodicity)','Color','w');
plot(x1_s, v1_s, 'k', 'LineWidth', 1.0); grid on;
xlabel('x_1'); ylabel('x_1'' (or v_1)');
title('Phase Portrait');

%% -----------------------------
% 9) 频谱：检查谐波含量（看 1× 是否主导）
% -----------------------------
figure('Name','Spectrum of steady-state','Color','w');
xw = x - mean(x);
N = numel(xw);
dt = mean(diff(t_s));
fs = 1/dt;

X = fft(xw .* hann(N));
f = (0:N-1)*(fs/N);
% 只画到 10×频率附近（足够看谐波）
f_max = min(max(f), 10*(Omega_test/(2*pi)));
Iplot = (f <= f_max);

plot(f(Iplot), 20*log10(abs(X(Iplot))/max(abs(X(Iplot))+eps) + eps), 'LineWidth', 1.2);
grid on;
xlabel('Frequency (Hz)'); ylabel('Normalized magnitude (dB)');
title('Spectrum (normalized) - check harmonics');

%% -----------------------------
% 10) （可选）如果你想同时把 FRF 曲线也画出来并标出 Omega_test
% -----------------------------
if FRF_ok
    figure('Name','FRF curve with Omega_test marker','Color','w');
    plot(Omega_frf, dB_frf, 'LineWidth', 1.5); grid on; hold on;
    xline(Omega_test, 'r--', 'LineWidth', 1.5);
    yline(dB_frf_1x, 'r:', 'LineWidth', 1.5);
    xlabel('\Omega'); ylabel('20log_{10}(|X_1|) (dB)');
    title('FRF (HBM) and verification point');
    legend('FRF','\Omega_{test}','FRF@ \Omega_{test}','Location','best');
end
