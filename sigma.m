%% 非线性验证 (倒序扫频 + 力传递率) - sigma 参数扫描
% 目标：倒序扫频 (Omega: 5.0 -> 0.2)，绘制力传递率 TF = |Ft|/Fw 的 dB 曲线
% sigma 取值为 [-1, 1, 1.2, 1.4, 1.6, 2.5]
% 横坐标：log scale (10 的幂次)

clear; clc; close all;

%% 基础参数设置（保持不变）
k1=1; k2=0.8;  %k1,k2可以调整  
L=4/9; U=2;    %非线性的几何参数
P.be1 = 1;     %上层刚度与（基准刚度也就是上层刚度之比）
P.mu  = 0.2;   %下层质量对上层质量的比值
P.be2 = 0.1;   %下层的等效刚度与基准刚度之比    P.be2 = 2*K2*(1-L)/L; 
P.al1 = -0.95; %上层准零刚度对应的等效线性刚度P.al1 = 2*K1*(1-L)/L; 
P.ga1 = k1/ (U^2 * L^3);   %上层准零刚度对应的三次方系数     
P.ga2 = k2/ (U^2 * L^3);   %下层准零刚度对应的三次方系数  
P.ze1 = 0.05;  %第二增阻尼比  需要根据质量比进行换算
% 电路参数（除sigma外）
P.lam   = 0.18;  %设置的第一层阻尼比为根据质量比进行计算，当电路断开时，模拟电磁分流阻尼比
P.kap_e = 1.0;   %电感系数
P.kap_c = 0.2;   %电容系数
% sigma_list = [-1, 1, 1.2, 1.4, 1.6, 2.5];  % 待扫描的电阻系数

% sigma_list = [-0.1, -0.5, -0.9, 1.2, 1.4, 1.6, 2.5];
sigma_list = [-0.1, -0.5, -0.9, 1.2, 1.4, 1.6, 2.5];  % 主要关注负电阻的隔振效果

%% 全局参数设置
global Fw FixedOmega
Fw = 0.005;         % 外激励幅值
FixedOmega = [];    % 扫频模式

%% 倒序扫频范围
Omega_Start = 5.0;
Omega_End   = 0.2;
Omega_Step  = -0.001;   % 负步长

fprintf('开始倒序扫频: %.2f -> %.2f\n', Omega_Start, Omega_End);
fprintf('sigma扫描值: ');
fprintf('%.1f ', sigma_list);
fprintf('\n\n');

%% 初始化图形
figure('Color','w', 'Position', [100, 100, 900, 600]);
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'XScale','log');
xlim(ax, [Omega_End Omega_Start]);
xlabel(ax, '\Omega (log scale)', 'FontSize', 12);
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)', 'FontSize', 12);
title(ax, '不同\sigma值下的力传递率曲线（负电阻增强隔振效果）', 'FontSize', 14);
% 0 dB 参考线
yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);

%% 颜色映射
colors = lines(length(sigma_list));

%% 对每个 sigma 进行扫频计算
for idx = 1:length(sigma_list)
    P.sigma = sigma_list(idx);  % 当前sigma值
    
    % 组装系统参数
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
    
    fprintf('正在计算 sigma = %.2f ... ', P.sigma);
    
    %% 高频起点求解（Newton）
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;  % 16维：[15维HB状态 + Omega]
    
    [x0_full, ok] = newton('nondim_temp2', y_init, sysP);
    if ~ok
        fprintf('高频起点求解失败！跳过该sigma值。\n');
        continue;
    end
    x0 = x0_full(1:15);
    
    %% 构造第二个点（启动弧长法）
    Omega_Next = Omega_Start + Omega_Step;
    y_init2 = [x0; Omega_Next];
    
    [x1_full, ok] = newton('nondim_temp2', y_init2, sysP);
    if ~ok
        fprintf('第二个点求解失败！跳过该sigma值。\n');
        continue;
    end
    x1 = x1_full(1:15);
    
    %% 弧长延拓（倒序扫频）
    [x_res, ~] = branch_follow2('nondim_temp2', 5000, Omega_Start, Omega_Next, x0, x1, sysP);
    
    %% 力传递率 TF（严格按文档 ft 定义：β2*x2 + γ2*x2^3 + 2µζ2*x2'）
    Om  = x_res(16,:).';     % Nx1
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);           % 代码里 ze1，但物理上就是 ζ2
    ga2 = sysP(11);          % γ2
    
    % x2 的HB系数：每一行 [x20, a21, b21, a23, b23]
    x2 = x_res(6:10,:).';    % Nx5
    
    % x2' 的HB系数（与 nondim_temp2 里 D(Ω) 一致）
    W = Om;
    x2_dot = zeros(size(x2));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2(:,3);
    x2_dot(:,3) = -W .* x2(:,2);
    x2_dot(:,4) = 3*W .* x2(:,5);
    x2_dot(:,5) = -3*W .* x2(:,4);
    
    % AFT：计算 x2^3 的HB系数（0/1/3）
    x2_cub = cubic_proj_013_batch(x2);   % Nx5
    
    % ft 的HB系数（0/1/3）
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;   % Nx5
    
    % 取基波幅值 + 三次谐波合成
    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);   % 1+3 合成
    
    TF    = ft_amp ./ Fw;
    TF_dB = 20*log10(max(TF, 1e-300));
    
    % 清理 NaN/Inf & 非正频率
    ok_idx = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_clean = Om(ok_idx);
    TF_dB_clean = TF_dB(ok_idx);
    
    % 绘制曲线
    if ~isempty(Om_clean)
        if P.sigma < 0
            % 负电阻用粗实线突出显示
            plot(ax, Om_clean, TF_dB_clean, '-', 'Color', colors(idx,:), ...
                 'LineWidth', 2.5, 'DisplayName', sprintf('\\sigma = %.2f (负电阻)', P.sigma));
        else
            % 正电阻用常规线宽
            plot(ax, Om_clean, TF_dB_clean, '-', 'Color', colors(idx,:), ...
                 'LineWidth', 1.8, 'DisplayName', sprintf('\\sigma = %.2f', P.sigma));
        end
        fprintf('完成，点数: %d\n', length(Om_clean));
    else
        fprintf('无有效数据\n');
    end
end

%% 图例和修饰
legend('Location', 'best', 'FontSize', 10);
grid on;
hold off;

%% ============ 本脚本内用的 AFT 函数（和 nondim_temp2 一致）===========
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