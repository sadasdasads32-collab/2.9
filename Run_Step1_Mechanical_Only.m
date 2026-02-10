%% Run_Step1_Mechanical_Only.m
% 审稿人特供：第一步 - 纯机械系统验证脚本
% 目的：
% 1. 令 lam=0，彻底切断电路对机械的反作用，验证纯机械双层隔振特性。
% 2. 检查 FRF.m 与 branch_follow2 的配合情况。

clear; clc; close all;

% --- 1. 定义系统参数 (SysP) ---
% 采用之前验证过的“强非线性”参数，但切断耦合
P.be1 = 1.0; 
P.be2 = 1.0; 
P.mu  = 0.5;
P.al1 = 0.0;   % 去掉线性QZS，突出非线性
P.ga1 = 0.5;   % 上层强非线性 (应当看到弯曲)
P.ze1 = 0.02;  % 低阻尼 (共振峰应当很尖)
P.lam = 0.0;   % <--- 【关键】设为0，断开机电耦合！
P.kap_e = 0.05; 
P.kap_c = 1.0; 
P.sigma = 0.5; 
P.ga2 = 0.2;   % 下层非线性

% 组装参数向量
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

% 全局变量设置
global Fw FixedOmega
Fw = 0.05;       % 激励力幅值
FixedOmega = []; % 告诉 nondim_temp2 这是扫频模式

fprintf('================================================\n');
fprintf('STEP 1: Mechanical Verification (Lambda = 0)\n');
fprintf('Target: Verify nonlinear peaks without electrical damping.\n');
fprintf('================================================\n');

% --- 2. 调用你的 FRF 扫频程序 ---
% 注意：FRF.m 内部硬编码了 Omega0=0.2, OmegaMax=6.0
% 如果需要修改范围，请直接去 FRF.m 里修改，或者将其改为参数输入
try
    % 调用 FRF.m
    % 返回的 x 是 16xN 矩阵，最后一行是频率 Omega
    x_results = FRF(sysP); 
catch ME
    fprintf('Error calling FRF: %s\n', ME.message);
    fprintf('请检查是否已将 newton.m, branch_follow2.m, nondim_temp2.m 放在同一目录下。\n');
    return;
end

% --- 3. 数据后处理与绘图 ---
if isempty(x_results)
    error('FRF returned empty results.');
end

Omega_vec = x_results(16, :); % 频率轴

% 提取幅值
Amp_X1 = sqrt(x_results(2,:).^2 + x_results(3,:).^2);
Amp_X2 = sqrt(x_results(7,:).^2 + x_results(8,:).^2);

% === 转 dB（相对 1）===
eps_db = 1e-12;                 % 防止 log(0)
Amp_X1_dB = 20*log10(Amp_X1 + eps_db);
Amp_X2_dB = 20*log10(Amp_X2 + eps_db);


% 绘图
figure('Position',[100,100,1000,400], 'Color', 'w');

subplot(1,2,1);
plot(Omega_vec, Amp_X1_dB, 'b-', 'LineWidth', 1.5);
xlabel('Frequency \Omega'); ylabel('Amplitude |X_1| (dB)');
title('Mechanical Response: Upper Mass (m1)'); grid on;

subplot(1,2,2);
plot(Omega_vec, Amp_X2_dB, 'r-', 'LineWidth', 1.5);
xlabel('Frequency \Omega'); ylabel('Amplitude |X_2| (dB)');
title('Mechanical Response: Lower Mass (m2)'); grid on;


fprintf('Done! Please inspect the curves.\n');
fprintf('Checklist:\n');
fprintf('1. Are there TWO peaks? (around 0.8 and 2.5)\n');
fprintf('2. Are they bending to the RIGHT? (Hardening)\n');
fprintf('3. Is the amplitude high? (Due to zero electrical damping)\n');