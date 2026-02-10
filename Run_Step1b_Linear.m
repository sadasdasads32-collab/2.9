%% Run_Step1b_Linear_Verification.m
% 审稿人特供：第二步验证 - 线性纯机械系统 (dB版 + 自动寻峰)
% 目的：
% 1. 验证线性模态 (Linear Modes)
% 2. 打印共振峰坐标，量化固有频率

clear; clc; close all;

% --- 1. 参数设置 (线性化) ---
P.be1 = 1.0;   
P.be2 = 1.0;   
P.mu  = 0.5;   
P.al1 = 0.0;   
P.ga1 = 0.0;   % 【关闭】上层非线性
P.ze1 = 0.05;  % 下层阻尼
P.lam = 0.0;   % 【断开】电路
P.kap_e = 0.05; 
P.kap_c = 1.0; 
P.sigma = 0.5; 
P.ga2 = 0.0;   % 【关闭】下层非线性

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.05;       
FixedOmega = []; 

fprintf('================================================\n');
fprintf('VERIFICATION 2: Linear Mechanical System (dB Mode)\n');
fprintf('================================================\n');

% --- 2. 调用 FRF ---
try
    x_results = FRF(sysP); 
catch ME
    fprintf('Error: %s\n', ME.message);
    return;
end

if isempty(x_results)
    error('No results.');
end

% --- 3. 数据处理 (转dB) ---
Omega_vec = x_results(16, :);

% 幅值 (Linear Scale)
Amp_X1_Lin = sqrt(x_results(2,:).^2 + x_results(3,:).^2);
Amp_X2_Lin = sqrt(x_results(7,:).^2 + x_results(8,:).^2);

% 幅值 (dB Scale)
Amp_X1_dB = 20 * log10(Amp_X1_Lin);
Amp_X2_dB = 20 * log10(Amp_X2_Lin);

% --- 4. 自动寻找共振峰 ---
[pks1, locs1] = find_peaks_custom(Amp_X1_dB, Omega_vec);
[pks2, locs2] = find_peaks_custom(Amp_X2_dB, Omega_vec);

fprintf('\n--- Resonance Peaks Report ---\n');
fprintf('X1 (Upper Mass) Peaks:\n');
for i = 1:length(pks1)
    fprintf('  Mode %d: Omega = %.4f, Amp = %.2f dB\n', i, locs1(i), pks1(i));
end

fprintf('X2 (Lower Mass) Peaks:\n');
for i = 1:length(pks2)
    fprintf('  Mode %d: Omega = %.4f, Amp = %.2f dB\n', i, locs2(i), pks2(i));
end
fprintf('------------------------------\n');

% --- 5. 绘图 (dB) ---
figure('Position',[100,100,1000,400], 'Color', 'w');

% X1
subplot(1,2,1);
plot(Omega_vec, Amp_X1_dB, 'b-', 'LineWidth', 2);
hold on;
% 标记峰值点
plot(locs1, pks1, 'ro', 'MarkerFaceColor', 'r'); 
xlabel('\Omega'); ylabel('|X_1| (dB)');
title('Linear Response: Upper Mass (m1)');
grid on; axis tight;

% X2
subplot(1,2,2);
plot(Omega_vec, Amp_X2_dB, 'r-', 'LineWidth', 2);
hold on;
plot(locs2, pks2, 'bo', 'MarkerFaceColor', 'b');
xlabel('\Omega'); ylabel('|X_2| (dB)');
title('Linear Response: Lower Mass (m2)');
grid on; axis tight;

fprintf('Done! Peaks are marked with circles.\n');

% --- 辅助函数：简单寻峰 (不依赖工具箱) ---
function [pks, locs] = find_peaks_custom(y, x)
    pks = [];
    locs = [];
    if length(y) < 3, return; end
    
    for i = 2:length(y)-1
        % 判断局部极大值
        if y(i) > y(i-1) && y(i) > y(i+1)
            % 简单的阈值过滤，防止噪音（可根据需要调整）
            if y(i) > -100 
                pks = [pks; y(i)];
                locs = [locs; x(i)];
            end
        end
    end
end