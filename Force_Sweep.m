%% Run_Auto_Force_Sweep.m
% 审稿人：Reviewer #2
% 功能：自动串联分析 —— 先扫频找峰，再定频扫力画 S 线
% 目的：证明系统在特定频率下的双稳态（S-curve）特性

clear; clc; close all;

% =========================================================================
% 1. 系统参数设置 (System Parameters)
% =========================================================================
% 选用强非线性参数以确保能看到 S 线
% sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
P.be1   = 1.0;   % sysP(1)
P.be2   = 1.0;   % sysP(2)
P.mu    = 0.2;   % sysP(3)
P.al1   = 1.5;   % sysP(4)
P.ga1   = 1.5;   % sysP(5)
P.ze1   = 0.05;  % sysP(6)
P.lam   = 0.05;

% 调谐后的电路参数（按你当前方案）
P.kap_e = 2.04;  % sysP(8)
P.kap_c = 1.0;   % sysP(9)
P.sigma = 0.5;   % sysP(10)
P.ga2   = 1.5;   % sysP(11)

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
        P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.005;          % 【FRF 激励幅值】：先用一个小一点的力扫频，或者中等力
FixedOmega = [];    % 确保 FRF 模式下为空

fprintf('================================================\n');
fprintf('STEP 1: Running Frequency Sweep (FRF) at Fw=%.4f\n', Fw);
fprintf('================================================\n');

% =========================================================================
% 2. 执行扫频 (Call FRF)
% =========================================================================
try
    % x_res 是 16xN 的矩阵，最后一行 x_res(16,:) 是频率 Omega
    x_res = FRF(sysP); 
catch ME
    error('FRF 调用失败，请检查 FRF.m 是否存在或报错: %s', ME.message);
end

if isempty(x_res)
    error('FRF 返回为空，无法进行后续分析。');
end

Omega_vec = x_res(16, :);
% 计算 X1 的幅值 (1次谐波)
Amp_X1 = sqrt(x_res(2,:).^2 + x_res(3,:).^2);

% 绘图：频响曲线
figure('Position', [100, 100, 1000, 400], 'Color', 'w');
subplot(1, 2, 1);
plot(Omega_vec, Amp_X1, 'k-', 'LineWidth', 1.2); hold on;
xlabel('\Omega'); ylabel('|X_1|'); title('Step 1: Frequency Response');
grid on;

% =========================================================================
% 3. 自动寻峰 (Find Peak)
% =========================================================================
[max_amp, idx_peak] = max(Amp_X1);
w_peak = Omega_vec(idx_peak);
x_state_peak = x_res(:, idx_peak); % 提取峰值处的完整状态向量 (16x1)

plot(w_peak, max_amp, 'ro', 'MarkerFaceColor', 'r');
text(w_peak, max_amp, sprintf('  Peak \\Omega=%.3f', w_peak), 'VerticalAlignment','bottom');
fprintf('   -> Found Resonance Peak at Omega = %.4f, Amp = %.4f\n', w_peak, max_amp);

% =========================================================================
% 4. 策略选择：定频扫力 (Force Sweep Strategy)
% =========================================================================
% 审稿人提示：
% 想要看到 S 线（多值性），通常需要在共振峰的“非线性影响区”扫力。
% 如果 ga1 > 0 (硬特性)，建议在 w_peak 或 w_peak 的右侧一点扫。
% 这里我们演示两个位置：
% Case A: 直接在峰值频率扫 (At Peak)
% Case B: 在峰值频率偏右 5% 处扫 (Post Peak) - 更有可能看到大幅度跳跃

target_Omega_list = [w_peak, w_peak * 1.1]; % 你可以修改这个系数
line_styles = {'b.-', 'm.-'};

subplot(1, 2, 2); hold on; grid on;
xlabel('Excitation Force F'); ylabel('|X_1|'); 
title('Step 2: Force Sweep (S-Curve Search)');

for k = 1:length(target_Omega_list)
    w_target = target_Omega_list(k);
    
    fprintf('\n------------------------------------------------\n');
    fprintf('STEP 2.%d: Sweeping Force at Omega = %.4f\n', k, w_target);
    fprintf('------------------------------------------------\n');
    
    % --- 构造 L1 的初值 ---
    % 这里的 trick 是：我们利用 FRF 在 w_peak 处的解 x_state_peak 作为初值猜测。
    % 如果 w_target 与 w_peak 相差不远，Newton迭代能拉过去。
    % 注意：L1 接受的 x_interest 是 [15x1 coeffs; Omega]
    
    x_guess = x_state_peak(1:15); % 仅取系数
    % 组合成 L1 需要的输入格式：[coeffs; w_target]
    x_input = [x_guess; w_target]; 
    
    % 设置扫力的初始力值。FRF 是在 Fw 下跑的，所以我们从 Fw 开始扫
    % L1(x_interest, sysP, StartForce/Step, BranchLen)
    % 第三个参数传 -0.001 表示步长（如果传值>0.1则被视为力），这里我们让它内部处理
    % 这里的 L1 逻辑是：如果 varargin{1} < F_switch，则视为步长。
    % 我们希望从 Fw 开始扫，Fw=0.02 可能小于阈值。
    % 让我们修改调用方式：
    % L1 内部: if abs(input) > F_switch(0.1) -> init force. 
    % Fw = 0.02 太小了，会被当成步长。
    % 既然我们要利用 x_input (它是 F=Fw 时的解)，L1 会自动把 x_input 的最后一行当作 Omega，
    % 但 L1 内部的 newton 需要一个力值。
    % 你的 L1 代码中: current_F = Fw (全局) 或 varargin{1}.
    % 我们显式设置全局 Fw 匹配当前状态
    Fw = 0.02; 
    
    % 调用 L1
    % 返回 branch_data: 16 x M，最后一行是 Force
    branch_data = L1(x_input, sysP, 0.001, 3000); 
    
    % --- 绘图 ---
    if ~isempty(branch_data)
        F_vals = branch_data(16, :);
        X1_vals = sqrt(branch_data(2,:).^2 + branch_data(3,:).^2);
        
        plot(F_vals, X1_vals, line_styles{k}, 'DisplayName', sprintf('\\Omega=%.3f', w_target));
        drawnow;
    else
        fprintf('   [Warning] L1 returned empty branch for Omega=%.3f\n', w_target);
    end
end

legend('show', 'Location', 'best');

fprintf('\nDone. Check the right subplot for S-curves.\n');