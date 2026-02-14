%% 非线性验证 - 电容参数扫描 (kap_c)
% 目标：倒序扫频 (Omega: 5.0 -> 0.2)，研究不同电容系数对力传递率的影响
% 电阻固定 sigma = 0.5（正电阻），只改变电容
% kap_c 扫描范围：[0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

clear; clc; close all;

%% 基础参数设置（保持不变）
k1=1; k2=0.8;  
L=4/9; U=2;    
P.be1 = 1;     
P.mu  = 0.2;   
P.be2 = 0.1;   
P.al1 = -0.95; 
P.ga1 = k1/ (U^2 * L^3);     
P.ga2 = k2/ (U^2 * L^3);     
P.ze1 = 0.05;  

% 电路参数
P.lam   = 0.18;      % 耦合强度
P.kap_e = 1.0;       % 电感系数（固定）
P.sigma = 0.5;       % ⭐ 电阻固定为 0.5（正电阻），不加负电阻
% P.sigma = -0.5;     % （备用选项，如果你想对比负电阻效果，取消注释这行）

% 电容扫描列表
kap_c_list = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0];

%% 全局参数设置
global Fw FixedOmega
Fw = 0.005;         
FixedOmega = [];    

%% 倒序扫频范围
Omega_Start = 5.0;
Omega_End   = 0.2;
Omega_Step  = -0.001;   

fprintf('========== 电容参数扫描 ==========\n');
fprintf('电阻固定: sigma = %.1f\n', P.sigma);
fprintf('电容扫描值: ');
fprintf('%.2f ', kap_c_list);
fprintf('\n\n');

%% 初始化图形 - 主图
figure('Color','w', 'Position', [100, 100, 1000, 700]);
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'XScale','log');
xlim(ax, [Omega_End Omega_Start]);   % 横轴固定扫频范围
% ylim(ax, [-60, 10]);   % 注释掉固定纵轴范围，让后续自动调整
xlabel(ax, '\Omega (log scale)', 'FontSize', 12);
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)', 'FontSize', 12);
title(ax, sprintf('不同电容系数 \\kappa_c 下的力传递率曲线 (\\sigma = %.1f, 正电阻)', P.sigma), 'FontSize', 14);
% 0 dB 参考线
yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

%% 颜色映射
colors = jet(length(kap_c_list));

%% 存储特定频率点的数据用于表格
freq_points = [1, 2, 3, 4, 5];  
table_data = zeros(length(kap_c_list), length(freq_points));

%% 对每个 kap_c 进行扫频计算
for idx = 1:length(kap_c_list)
    P.kap_c = kap_c_list(idx);  % 当前电容值
    
    % 组装系统参数
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
    
    fprintf('正在计算 kap_c = %.3f ... ', P.kap_c);
    
    %% 高频起点求解
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;  
    
    [x0_full, ok] = newton('nondim_temp2', y_init, sysP);
    if ~ok
        fprintf('高频起点求解失败！跳过该kap_c值。\n');
        continue;
    end
    x0 = x0_full(1:15);
    
    %% 构造第二个点
    Omega_Next = Omega_Start + Omega_Step;
    y_init2 = [x0; Omega_Next];
    
    [x1_full, ok] = newton('nondim_temp2', y_init2, sysP);
    if ~ok
        fprintf('第二个点求解失败！跳过该kap_c值。\n');
        continue;
    end
    x1 = x1_full(1:15);
    
    %% 弧长延拓
    [x_res, ~] = branch_follow2('nondim_temp2', 5000, Omega_Start, Omega_Next, x0, x1, sysP);
    
    %% 力传递率计算
    Om  = x_res(16,:).';     
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);           
    ga2 = sysP(11);          
    
    % x2 的HB系数
    x2 = x_res(6:10,:).';    
    
    % x2' 的HB系数
    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);
    
    % AFT：计算 x2^3
    x2_cub = cubic_proj_013_batch(x2);   
    
    % ft 的HB系数
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;   
    
    % 合成幅值
    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);   
    
    TF    = ft_amp ./ Fw;
    TF_dB = 20*log10(max(TF, 1e-300));
    
    % 清理无效数据
    ok_idx = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_clean = Om(ok_idx);
    TF_dB_clean = TF_dB(ok_idx);
    
    % 绘制曲线
    if ~isempty(Om_clean)
        plot(ax, Om_clean, TF_dB_clean, '-', 'Color', colors(idx,:), ...
             'LineWidth', 2.0, 'DisplayName', sprintf('\\kappa_c = %.3f', P.kap_c));
        
        % 提取特定频率点的数据
        for j = 1:length(freq_points)
            [~, nearest_idx] = min(abs(Om_clean - freq_points(j)));
            if ~isempty(nearest_idx)
                table_data(idx, j) = TF_dB_clean(nearest_idx);
            end
        end
        
        fprintf('完成，点数: %d\n', length(Om_clean));
    else
        fprintf('无有效数据\n');
    end
end

%% 自动调整纵轴范围以完整显示所有曲线
ax = gca;
axis tight;
y_limits = ylim(ax);
y_margin = 0.05 * (y_limits(2) - y_limits(1));
ylim(ax, [y_limits(1) - y_margin, y_limits(2) + y_margin]);

%% 图例和修饰
legend('Location', 'eastoutside', 'FontSize', 9);
grid on;
hold off;
%% 绘制第二个图：特定频率点的电容影响曲线
figure('Color','w', 'Position', [150, 150, 900, 500]);
ax2 = gca; hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

% 绘制不同频率下TF_dB随kap_c的变化
for j = 1:length(freq_points)
    plot(ax2, kap_c_list, table_data(:,j), 'o-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', sprintf('\\Omega = %d', freq_points(j)));
end

xlabel(ax2, '\kappa_c (电容系数)', 'FontSize', 12);
ylabel(ax2, 'Force Transmissibility (dB)', 'FontSize', 12);
title(ax2, sprintf('不同频率下力传递率随电容系数的变化 (\\sigma = %.1f, 正电阻)', P.sigma), 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
set(ax2, 'XScale', 'log');  % 电容取对数坐标
xlim(ax2, [min(kap_c_list)*0.9, max(kap_c_list)*1.1]);
grid on;

%% 输出表格数据
fprintf('\n\n===== 力传递率数据表 (dB) =====\n');
fprintf('电阻固定: σ = %.1f\n', P.sigma);
fprintf('Ω\t');
for j = 1:length(freq_points)
    fprintf('%.1f\t', freq_points(j));
end
fprintf('\n');
fprintf('----------------------------------------\n');
for idx = 1:length(kap_c_list)
    fprintf('κ_c=%.2f\t', kap_c_list(idx));
    for j = 1:length(freq_points)
        if table_data(idx, j) ~= 0
            fprintf('%.1f\t', table_data(idx, j));
        else
            fprintf('N/A\t');
        end
    end
    fprintf('\n');
end

%% 如果想对比负电阻效果，取消下面的注释块
%{
%% 对比：负电阻情况 (sigma = -0.5)
fprintf('\n\n========== 负电阻对比 ==========\n');
P.sigma = -0.5;  % 临时改为负电阻
... 重复上述计算 ...
%}

%% ============ AFT 函数 ===========
function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    X_time  = (T_mat * U.').';       
    X3_time = X_time.^3;
    cubic   = (T_inv * X3_time.').'; 
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