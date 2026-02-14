%% 纯机械双层QZS系统时域验证 (增强版)
clear; clc; close all;

%% 1. 系统参数
k1 = 1; k2 = 0.8;
L = 4/9; U = 2;
P.be1 = 1;
P.mu  = 0.2;
P.be2 = 0.1;
P.al1 = -0.95;
P.ga1 = k1 / (U^2 * L^3);
P.ga2 = k2 / (U^2 * L^3);
P.ze1 = 0.05;

% 电路全零
P.lam   = 0.0;
P.kap_e = 0.0;
P.kap_c = 0.0;
P.sigma = 0.0;

Fw = 0.005;
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

%% 2. 待验证频率点
Omega_list = [1, 2, 3, 4, 5];  % 避免Ω=0
n_freq = length(Omega_list);

% 预分配
TF_fd = zeros(1, n_freq);
TF_td = zeros(1, n_freq);
x_fd_all = cell(1, n_freq);  % 存储频域解，作为时域初值

%% 3. 频域逐点求解
fprintf('========== 频域逐点 Newton ==========\n');
global Fw FixedOmega
Fw = Fw;
FixedOmega = [];

for idx = 1:n_freq
    Omega = Omega_list(idx);
    fprintf('Ω = %.2f ... ', Omega);
    
    % 使用更好的初值：从低频解延续（对Ω=1用扰动，后续用前一个解）
    if idx == 1
        y0 = zeros(15,1);
        y0(2) = 0.1;  % 稍大的扰动，帮助起振
    else
        y0 = x_fd_all{idx-1}(1:15);  % 使用上一个频率的解作为初值
    end
    y0_full = [y0; Omega];
    
    [x_full, ok, ~] = newton('nondim_temp2', y0_full, sysP);
    
    if ~ok
        fprintf('Newton失败\n');
        TF_fd(idx) = NaN;
        continue;
    end
    
    x_fd_all{idx} = x_full;
    
    % 计算力传递率
    x2 = x_full(6:10);
    be2 = P.be2; ga2 = P.ga2; mu = P.mu; ze2 = P.ze1;
    W = Omega;
    x2_dot = [0; W*x2(3); -W*x2(2); 3*W*x2(5); -3*W*x2(4)];
    x2_cub = cubic_proj_013(x2);
    ft_coef = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
    ft1 = sqrt(ft_coef(2)^2 + ft_coef(3)^2);
    ft3 = sqrt(ft_coef(4)^2 + ft_coef(5)^2);
    ft_amp = sqrt(ft1^2 + ft3^2);
    TF_fd(idx) = ft_amp / Fw;
    
    fprintf('TF = %.2f dB\n', 20*log10(TF_fd(idx)));
end

%% 4. 时域仿真参数
t_total    = 5000;      % 总时间
n_cycles   = 50;        % 用于分析的周期数
fprintf('\n========== 时域仿真 (ode45) ==========\n');

for idx = 1:n_freq
    Omega = Omega_list(idx);
    fprintf('Ω = %.2f ... ', Omega);
    
    % 从频域解构造时域初值（关键改进！）
    if ~isempty(x_fd_all{idx})
        x_full = x_fd_all{idx};
        % 从HB系数重构时域初始状态
        [~, T_mat, ~] = get_AFT_matrices();
        
        % ξ1 及其导数
        x1_hb = x_full(1:5);
        x1_time = T_mat * x1_hb;
        x1_dot_hb = [0; Omega*x1_hb(3); -Omega*x1_hb(2); 3*Omega*x1_hb(5); -3*Omega*x1_hb(4)];
        x1_dot_time = T_mat * x1_dot_hb;
        
        % ξ2 及其导数
        x2_hb = x_full(6:10);
        x2_time = T_mat * x2_hb;
        x2_dot_hb = [0; Omega*x2_hb(3); -Omega*x2_hb(2); 3*Omega*x2_hb(5); -3*Omega*x2_hb(4)];
        x2_dot_time = T_mat * x2_dot_hb;
        
        % Q 及其导数（电路全零，设为0）
        Q_time = zeros(size(x1_time));
        Q_dot_time = zeros(size(x1_time));
        
        % 取t=0时刻的值作为初值
        y0 = [x1_time(1); x1_dot_time(1); x2_time(1); x2_dot_time(1); Q_time(1); Q_dot_time(1)];
    else
        y0 = zeros(6,1);  % 如果没有频域解，用零初值
        y0(2) = 0.1;      % 加扰动
    end
    
    % 定义odefun
    odefun = @(tau, y) odeFunc(tau, y, P, Fw, Omega);
    
    % 求解
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    sol = ode45(odefun, [0, t_total], y0, opts);
    
    % 提取稳态段
    T_period = 2*pi/Omega;
    t_start = t_total - n_cycles * T_period;
    t_eval = linspace(t_start, t_total, 10000);  % 更多采样点
    y_eval = deval(sol, t_eval);
    
    % 绘制时域响应检查（可选）
    if idx == 1  % 只对第一个频率绘图检查
        figure(100);
        subplot(2,1,1);
        plot(t_eval, y_eval(3,:));  % ξ2
        xlabel('\tau'); ylabel('\xi_2');
        title(sprintf('Ω = %.2f 时域响应', Omega));
        grid on;
        subplot(2,1,2);
        plot(t_eval(end-500:end), y_eval(3,end-500:end));  % 最后部分
        xlabel('\tau'); ylabel('\xi_2 (稳态)');
        grid on;
    end
    
    xi2 = y_eval(3,:);
    xi2_dot = y_eval(4,:);
    
    % 计算ft
    ft = P.be2 * xi2 + P.ga2 * xi2.^3 + 2 * P.mu * P.ze1 * xi2_dot;
    
    % FFT分析
    N = length(ft);
    dt = t_eval(2) - t_eval(1);
    ft_fft = fft(ft);
    freq = (0:N-1) / (N*dt);
    
    % 找到精确频率点（使用插值提高精度）
    [~, i1] = min(abs(freq - Omega));
    % 使用附近点加权平均
    if i1 > 1 && i1 < N
        w = [0.5, 1, 0.5];
        idx1 = i1-1:i1+1;
        amp1 = sqrt(sum(w .* (abs(ft_fft(idx1)).^2)) / sum(w)) * 2 / N;
    else
        amp1 = 2 * abs(ft_fft(i1)) / N;
    end
    
    [~, i3] = min(abs(freq - 3*Omega));
    if i3 > 1 && i3 < N
        w = [0.5, 1, 0.5];
        idx3 = i3-1:i3+1;
        amp3 = sqrt(sum(w .* (abs(ft_fft(idx3)).^2)) / sum(w)) * 2 / N;
    else
        amp3 = 2 * abs(ft_fft(i3)) / N;
    end
    
    amp_total = sqrt(amp1^2 + amp3^2);
    TF_td(idx) = amp_total / Fw;
    
    fprintf('TF = %.2f dB\n', 20*log10(TF_td(idx)));
end

%% 5. 绘制对比图
figure('Color','w','Position',[100,100,900,600]);
semilogx(Omega_list, 20*log10(abs(TF_fd)), 'bo-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(Omega_list, 20*log10(abs(TF_td)), 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('\Omega (log scale)','FontSize',12);
ylabel('Force Transmissibility (dB)','FontSize',12);
title('纯机械双层QZS系统：频域 vs 时域 (初值取自频域解)','FontSize',14);
legend('频域 (逐点Newton)','时域 (ode45)','Location','best');
grid on; 
xlim([min(Omega_list)*0.9, max(Omega_list)*1.1]);
set(gca,'XScale','log');

%% 6. 输出表格
fprintf('\n\n===== 力传递率对比表 (dB) =====\n');
fprintf('Ω\t频域\t时域\t差值\n');
for idx = 1:n_freq
    fprintf('%.1f\t%.2f\t%.2f\t%.2f\n', ...
        Omega_list(idx), 20*log10(TF_fd(idx)), 20*log10(TF_td(idx)), ...
        20*log10(TF_fd(idx)) - 20*log10(TF_td(idx)));
end

%% ============ 辅助函数 ============
function dydt = odeFunc(tau, y, P, Fw, Omega)
    xi1 = y(1); xi1_dot = y(2);
    xi2 = y(3); xi2_dot = y(4);
    Q   = y(5); Q_dot   = y(6);
    
    be1 = P.be1; al1 = P.al1; ga1 = P.ga1;
    be2 = P.be2; ga2 = P.ga2;
    mu = P.mu; ze2 = P.ze1;
    
    x12 = xi1 - xi2;
    upper_force = (be1+al1)*x12 + ga1*x12^3;
    F_ext = Fw * cos(Omega * tau);
    
    % 纯机械（电路全零）
    xi1_ddot = F_ext - upper_force;
    xi2_ddot = ( -2*mu*ze2*xi2_dot - be2*xi2 - ga2*xi2^3 + upper_force ) / mu;
    
    dydt = [xi1_dot; xi1_ddot; xi2_dot; xi2_ddot; 0; 0];
end

function cubic = cubic_proj_013(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    cubic = T_inv * ((T_mat * u).^3);
end

function [N, T_mat, T_inv] = get_AFT_matrices()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64;
        t = (0:pN-1)' * (2*pi/pN);
        c1 = cos(t); s1 = sin(t);
        c3 = cos(3*t); s3 = sin(3*t);
        dc = ones(pN,1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end
    N = pN; T_mat = pT; T_inv = pTinv;
end