%% 带电路单组参数验证 - 弧长延拓修正版（双曲线对比） 含性能指标对比
% 依赖外部函数: nondim_temp2, newton, branch_follow2

clc; clear; close all;

%% -------- 1. 基础参数定义（机械部分）--------
% Wang 图注给定参数（BG Model）
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.0;       % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 反推 v 与非线性系数
v = 2.5;        % 由 L=4/9, K1=1 反推
alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

% 机械部分固定参数
P_base.be1 = 1.0;
P_base.al1 = alpha1 - P_base.be1;
P_base.be2 = alpha2;
P_base.ga1 = gamma1;
P_base.ga2 = gamma2;
P_base.mu  = mu;
P_base.ze1 = 0.05;   % 下层阻尼比

global Fw
Fw = 0.005;

%% -------- 2. 定义两组电路参数（待对比）--------
% 原参数组
P1 = P_base;
P1.lam   = 0.18;
P1.kap_e = 0.395;
P1.kap_c = 0.032;
P1.sigma = 0.623;

% 新参数组（纯机械，无电路影响）
P2 = P_base;
P2.lam   = 0;
P2.kap_e = 0;
P2.kap_c = 0;
P2.sigma = 0;

% 将两组参数放入元胞数组，便于循环处理
param_list = {P1, P2};
names = {'NIC电路 (\lambda=0.18, \kappa_e=1.82, \kappa_c=0.01, \sigma=0.54)', ...
         '最优参数对比曲线 (\lambda=0, \kappa_e=0, \kappa_c=0, \sigma=0)'};
colors = {'b', 'r'};   % 蓝色和红色

%% -------- 3. 弧长延拓通用设置 --------
Omega_Start = 10.0;
Omega_Step  = -0.01;  % 初始步长
Omega_Next  = Omega_Start + Omega_Step;

% 初始化图形
figure('Color','w', 'Position',[150 150 700 500]);
ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log');
xlabel(ax, '\Omega (log scale)');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title(ax, 'BG Model 电路参数对比');
yline(ax, 0, 'k--', '0 dB');

%% -------- 4. 准备存储性能指标 --------
% 使用结构体数组保存每条曲线的结果
results = struct('name', {}, 'peak_TF_dB', {}, 'peak_disp', {}, ...
                 'f_cross_0dB', {}, 'f_cross_m40dB', {});

%% -------- 5. 循环计算并绘制两条曲线 --------
for k = 1:2
    P = param_list{k};
    
    % 组装系统参数向量 sysP
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
    
    fprintf('\n========================================\n');
    fprintf('开始计算第 %d 组: %s\n', k, names{k});
    fprintf('参数: lam=%.2f, kap_e=%.2f, kap_c=%.2f, sigma=%.2f\n', ...
            P.lam, P.kap_e, P.kap_c, P.sigma);
    
    %% -------- 5.1 求解初始两个点 --------
    % 第一个点 (高频起点)
    y_init = zeros(15,1);
    y_init(end+1) = Omega_Start;
    [x0_full, ok0] = newton('nondim_temp2', y_init, sysP);
    
    if ~ok0
        warning('第 %d 组: 高频起点求解失败，跳过该曲线。', k);
        continue;
    end
    x0 = x0_full(1:15);
    
    % 第二个点 (用于确定弧长初始切线方向)
    y_init2 = [x0; Omega_Next];
    [x1_full, ok1] = newton('nondim_temp2', y_init2, sysP);
    
    if ~ok1
        warning('第 %d 组: 第二个初始点求解失败，跳过该曲线。', k);
        continue;
    end
    x1 = x1_full(1:15);
    
    %% -------- 5.2 弧长延拓主循环 --------
    [x_res, ~] = branch_follow2('nondim_temp2', 3000, Omega_Start, Omega_Next, x0, x1, sysP);
    
    %% -------- 5.3 提取响应数据 --------
    Om  = x_res(16,:).';                % 频率列向量
    x1_h = x_res(1:5,:).';              % 上层位移谐波系数 [dc, cos1, sin1, cos3, sin3]
    x2_h = x_res(6:10,:).';             % 下层位移谐波系数 [dc, cos1, sin1, cos3, sin3]
    
    be2 = sysP(2);
    mu  = sysP(3);
    ze2 = sysP(6);
    ga2 = sysP(11);
    
    % 计算 x2_dot
    W = Om;
    x2_dot = zeros(size(x2_h));
    x2_dot(:,1) = 0;
    x2_dot(:,2) = W .* x2_h(:,3);
    x2_dot(:,3) = -W .* x2_h(:,2);
    x2_dot(:,4) = 3*W .* x2_h(:,5);
    x2_dot(:,5) = -3*W .* x2_h(:,4);
    
    % 批量计算三次非线性项投影
    x2_cub = cubic_proj_013_batch(x2_h);
    
    % 合成传递力
    ft = be2*x2_h + ga2*x2_cub + 2*mu*ze2*x2_dot;
    ft1 = hypot(ft(:,2), ft(:,3));
    ft3 = hypot(ft(:,4), ft(:,5));
    ft_amp = hypot(ft1, ft3);
    
    % 转换为分贝
    TF    = ft_amp ./ Fw;
    TF_dB = 20*log10(max(TF, 1e-300));
    
    % 剔除无效点
    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_valid = Om(valid);
    TF_dB_valid = TF_dB(valid);
    
    %% -------- 5.4 绘制当前曲线 --------
    plot(ax, Om_valid, TF_dB_valid, 'Color', colors{k}, 'LineWidth', 1.5, ...
         'DisplayName', names{k});
    
    %% -------- 5.5 计算性能指标 --------
    % 上层位移基频幅值
    amp1 = sqrt(x1_h(:,2).^2 + x1_h(:,3).^2);   % 忽略直流和高次谐波
    % 注意：有效点筛选应与TF一致
    amp1 = amp1(valid);
    
    % 排序（频率升序）
    [Om_sorted, sort_idx] = sort(Om_valid);
    TF_sorted = TF_dB_valid(sort_idx);
    amp1_sorted = amp1(sort_idx);
    
    % 峰值力传递率
    peak_TF = max(TF_sorted);
    
    % 峰值动态位移（上层）
    peak_disp = max(amp1_sorted);
    
    % 0dB穿越频率
    f_cross_0 = find_cross_frequency(Om_sorted, TF_sorted, 0);
    
    % -40dB穿越频率
    f_cross_m40 = find_cross_frequency(Om_sorted, TF_sorted, -40);
    
    %% -------- 5.6 存储结果 --------
    results(k).name = names{k};
    results(k).peak_TF_dB = peak_TF;
    results(k).peak_disp = peak_disp;
    results(k).f_cross_0dB = f_cross_0;
    results(k).f_cross_m40dB = f_cross_m40;
    
    fprintf('第 %d 组计算完成。\n', k);
end

%% -------- 6. 输出性能指标对比表格 --------
fprintf('\n\n==================== 性能指标对比 ====================\n');
fprintf('%-55s %-16s %-16s %-16s %-16s\n', ...
        '曲线名称', '峰值力传递率(dB)', '峰值位移', '0dB穿越频率', '-40dB穿越频率');
fprintf('----------------------------------------------------------------------------------------\n');
for k = 1:2
    if isempty(results(k).name)
        continue;
    end
    name = results(k).name;
    peak_TF = results(k).peak_TF_dB;
    peak_disp = results(k).peak_disp;
    f0 = results(k).f_cross_0dB;
    fm40 = results(k).f_cross_m40dB;
    
    % 处理可能的NaN（无穿越）
    if isnan(f0)
        f0_str = '无';
    else
        f0_str = sprintf('%.4f', f0);
    end
    if isnan(fm40)
        fm40_str = '无';
    else
        fm40_str = sprintf('%.4f', fm40);
    end
    
    fprintf('%-55s %-16.4f %-16.4f %-16s %-16s\n', ...
            name, peak_TF, peak_disp, f0_str, fm40_str);
end
fprintf('========================================================\n');

%% -------- 7. 图形修饰 --------
legend(ax, 'Location', 'best');
xlim(ax, [0.1, Omega_Start]);
hold(ax,'off');

%% ============ 辅助函数 ============
function f_cross = find_cross_frequency(f, y, level)
    % 寻找第一次从上方穿越给定水平线的频率
    % 输入：
    %   f : 频率升序数组
    %   y : 对应的传递率 dB 值
    %   level : 目标 dB 值
    % 输出：
    %   f_cross : 穿越频率，若无穿越则返回 NaN

    % 找到第一个 <= level 的点
    idx = find(y <= level, 1, 'first');
    if isempty(idx)
        f_cross = NaN;
        return;
    end
    
    if idx == 1
        % 第一个点就已经低于 level，穿越发生在最小频率之前
        % 保守估计取最小频率作为穿越频率
        f_cross = f(1);
    else
        % 在前一点和当前点之间线性插值
        f1 = f(idx-1); y1 = y(idx-1);
        f2 = f(idx);   y2 = y(idx);
        if y1 == y2
            f_cross = f1;
        else
            f_cross = f1 + (level - y1) * (f2 - f1) / (y2 - y1);
        end
    end
end

function cubic = cubic_proj_013_batch(U)
    % U : N x 5 矩阵，每行为一组谐波系数 [dc, cos1, sin1, cos3, sin3]
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    X_time  = (T_mat * U.').';       % 时域转换
    X3_time = X_time.^3;             % 立方非线性
    cubic   = (T_inv * X3_time.').'; % 频域投影
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