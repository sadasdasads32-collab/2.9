%% Run_Optimization_Map.m
% 双参数优化扫描：寻找最佳 RLC 参数 (kap_e, sigma)
% 目标：最小化最大共振峰值 (Minimax Optimization)
% 
% X轴：电感 kap_e (影响调谐频率)
% Y轴：电阻 sigma (影响电磁阻尼)
% Z轴：Max(|X1|) dB

clear; clc; close all;

%% 1. 设置扫描范围 (根据你的系统特性调整)
% 建议先粗扫，锁定范围后再细扫
kape_vec = linspace(0.01, 3.0, 25);  % 电感范围
sigma_vec = linspace(0.1, 3.0, 25);  % 电阻范围

[K_grid, S_grid] = meshgrid(kape_vec, sigma_vec);
Z_peak_dB = zeros(size(K_grid));     % 存储最大峰值

%% 2. 固定参数设置
global Fw FixedOmega
Fw = 0.05; 
FixedOmega = [];

% 基础参数 (保持 Lambda = 0.5 不变，专门优化 RLC)
P_base.be1 = 1.0; P_base.be2 = 1.0; P_base.mu = 0.2;
P_base.al1 = 0.0; P_base.ga1 = 0.1; P_base.ze1 = 0.05;
P_base.lam = 0.5; % <--- 固定在强耦合工况
P_base.kap_c = 1.0; 
P_base.ga2 = 0.2; 

%% 3. 开始双重循环扫描
total_steps = numel(K_grid);
fprintf('开始参数优化扫描，共 %d 个工况...\n', total_steps);

tic;
for i = 1:size(K_grid, 1)     % sigma loop
    for j = 1:size(K_grid, 2) % kap_e loop
        
        % 当前参数
        curr_kap_e = K_grid(i,j);
        curr_sigma = S_grid(i,j);
        
        % 组装 sysP
        sysP = [P_base.be1; P_base.be2; P_base.mu; P_base.al1; ...
                P_base.ga1; P_base.ze1; P_base.lam; ...
                curr_kap_e; P_base.kap_c; curr_sigma; P_base.ga2];
        
        % --- 核心：调用 FRF 获取曲线 ---
        try
            % 这里的 FRF 需要能够静默运行（不画图）
            % 如果 FRF 内部有 plot 命令，建议注释掉或加个开关
            x_res = FRF(sysP); 
            
            if isempty(x_res)
                Z_peak_dB(i,j) = NaN; % 计算失败
            else
                % 提取幅值曲线
                Omega = x_res(16, :);
                AmpX1 = sqrt(x_res(2,:).^2 + x_res(3,:).^2);
                
                % 提取最大峰值 (H_inf norm)
                max_amp = max(AmpX1);
                Z_peak_dB(i,j) = 20*log10(max_amp + eps);
                
                % (进阶) 可以在这里加稳定性判断：
                % 如果曲线出现极大的跳跃或不连续，视为不稳定，设为 NaN
            end
        catch
            Z_peak_dB(i,j) = NaN;
        end
        
        % 进度条
        current_idx = (i-1)*size(K_grid,2) + j;
        if mod(current_idx, 10) == 0
            fprintf('进度: %.1f%% (Peak = %.2f dB)\n', current_idx/total_steps*100, Z_peak_dB(i,j));
        end
    end
end
toc;

%% 4. 绘图：优化云图 (Contour Plot)
figure('Name','Optimization Map','Color','w','Position',[100,100,600,500]);

% 绘制等高线填充图
% 颜色越蓝 -> 峰值越小 -> 效果越好
contourf(K_grid, S_grid, Z_peak_dB, 20, 'LineStyle','none'); 
hold on;
colormap(jet); % 或者 parula
colorbar;
caxis([-10, 5]); % 根据你的数据范围手动调整颜色映射范围，突出最小值

% 标出全局最优点
[min_val, min_idx] = min(Z_peak_dB(:));
[r_opt, c_opt] = ind2sub(size(Z_peak_dB), min_idx);
opt_kap = K_grid(r_opt, c_opt);
opt_sig = S_grid(r_opt, c_opt);

plot(opt_kap, opt_sig, 'kp', 'MarkerSize', 15, 'MarkerFaceColor','w');
text(opt_kap, opt_sig, sprintf('  Optimal\n  %.2f dB', min_val), 'Color','w','FontWeight','bold');

% 装饰
xlabel('Inductance Parameter \kappa_e (kap\_e)');
ylabel('Resistance Parameter \sigma (sigma)');
title(['H_{\infty} Optimization Map (\lambda = ' num2str(P_base.lam) ')']);
subtitle('Color represents Max Resonance Peak (dB)');

% --- 你的“研究发现”将在这里体现 ---
% 观察：
% 1. 是否有一个明显的“深蓝山谷”？
% 2. 最优点是否对应“两个峰等高”的情况？