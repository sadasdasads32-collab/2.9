%% 多目标电路参数优化脚本：基于 Pareto 遗传算法 (gamultiobj)
% 目标 1: 最小化峰值力传递率 (dB)
% 目标 2: 最小化峰值动态位移 (m1)
% 目标 3: 最大化单位隔离频带 (即最小化 0dB 穿越频率)
% 目标 4: 最大化 -40dB 隔离频带 (即最小化 -40dB 穿越频率)
%
% 依赖外部函数: nondim_temp2, newton, branch_follow2

clc; clear; close all;

%% -------- 1. 基础参数与固定参数定义 --------
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0;       % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 反推 v 与非线性系数
v = 2.5;        
alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

P.be1 = 1.0;
P.al1 = alpha1 - P.be1;
P.be2 = alpha2;
P.ga1 = gamma1;
P.ga2 = gamma2;
P.mu  = mu;
P.ze1 = 0.05;   

P.lam = 0.18;   % 耦合系数固定

global Fw
Fw = 0.008;

% 组装系统参数向量模板 [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
sysP_template = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, 0, 0, 0, P.ga2];

Omega_Start = 10.0;
Omega_Step  = -0.01; 
Omega_Next  = Omega_Start + Omega_Step;
Nsteps      = 3000;

%% -------- 2. 优化设置 (帕累托多目标) --------
% 变量 x = [sigma, kap_e, kap_c]
lb = [0.01, 0.01, 0.01];      % 下界 (避免纯0导致的刚性矩阵奇异)
ub = [2.0,  2.0,  2.0 ];      % 上界

% 帕累托算法配置
options = optimoptions('gamultiobj', ...
    'Display', 'iter', ...           % 打印迭代
    'PopulationSize', 60, ...        % 种群大小 (多目标建议设大些)
    'MaxGenerations', 40, ...        % 最大进化代数
    'UseParallel', true, ...         % 开启并行池加快计算
    'ParetoFraction', 0.35, ...      % 帕累托前沿保留比例
    'FunctionTolerance', 1e-4);

obj_fun = @(x) multiobj_func(x, sysP_template, Omega_Start, Omega_Next, Nsteps);

%% -------- 3. 运行多目标遗传算法 --------
fprintf('启动 NSGA-II 帕累托多目标优化...\n');
fprintf('优化的 4 个目标: [峰值传递率(dB), 峰值位移(m1), 0dB穿越频率, -40dB穿越频率]\n');
tic;
[x_pareto, fval_pareto] = gamultiobj(obj_fun, 3, [], [], [], [], lb, ub, [], options);
toc;

%% -------- 4. 结果后处理：提取折中解 (Utopia Point Method) --------
% 帕累托前沿是一组解，我们需要自动挑出一个 4 个目标都比较均衡的“最优折中解”
% 对适应度进行归一化
f_min = min(fval_pareto, [], 1);
f_max = max(fval_pareto, [], 1);
f_norm = (fval_pareto - f_min) ./ max(f_max - f_min, 1e-6);

% 计算每个解到“理想乌托邦点 (0,0,0,0)”的欧氏距离
distances = sum(f_norm.^2, 2);
[~, best_idx] = min(distances);

x_opt = x_pareto(best_idx, :);
f_opt = fval_pareto(best_idx, :);

fprintf('\n========== 帕累托折中解 (Utopia Compromise) ==========\n');
fprintf('最优电路参数:\n');
fprintf('  sigma (电阻) = %.4f\n', x_opt(1));
fprintf('  kap_e (电感) = %.4f\n', x_opt(2));
fprintf('  kap_c (电容) = %.4f\n', x_opt(3));
fprintf('对应的 4 大性能指标:\n');
fprintf('  1. 峰值力传递率 = %.2f dB\n', f_opt(1));
fprintf('  2. 峰值动态位移 = %.4f\n', f_opt(2));
fprintf('  3. 0dB 穿越频率 = %.2f (越小频带越宽)\n', f_opt(3));
fprintf('  4. -40dB 穿越频率 = %.2f (越小频带越宽)\n', f_opt(4));

%% -------- 5. 绘制折中解的验证频响曲线 --------
sysP_opt = sysP_template;
sysP_opt(10) = x_opt(1); % sigma
sysP_opt(8)  = x_opt(2); % kap_e
sysP_opt(9)  = x_opt(3); % kap_c

y_init = [zeros(15,1); Omega_Start];
[x0_full, ~] = newton('nondim_temp2', y_init, sysP_opt);
x0 = x0_full(1:15);
y_init2 = [x0; Omega_Next];
[x1_full, ~] = newton('nondim_temp2', y_init2, sysP_opt);
x1 = x1_full(1:15);
[x_res, ~] = branch_follow2('nondim_temp2', Nsteps, Omega_Start, Omega_Next, x0, x1, sysP_opt);

Om  = x_res(16,:).';
be2 = sysP_opt(2); mu = sysP_opt(3); ze2 = sysP_opt(6); ga2 = sysP_opt(11);
x2 = x_res(6:10,:).';

W = Om;
x2_dot = zeros(size(x2));
x2_dot(:,2) = W .* x2(:,3); x2_dot(:,3) = -W .* x2(:,2);
x2_dot(:,4) = 3*W .* x2(:,5); x2_dot(:,5) = -3*W .* x2(:,4);
x2_cub = cubic_proj_013_batch(x2);

ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
ft1 = hypot(ft(:,2), ft(:,3)); ft3 = hypot(ft(:,4), ft(:,5));
ft_amp = hypot(ft1, ft3);
TF_dB = 20*log10(max(ft_amp ./ Fw, 1e-300));

valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);

figure('Color','w','Position',[150 150 700 500]);
ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log');
xlabel(ax,'\Omega (log scale)');
ylabel(ax,'Force Transmissibility (dB)');
title(sprintf('Pareto Optimal: \\sigma=%.2f, \\kappa_e=%.2f, \\kappa_c=%.2f\nTF_{peak}=%.1fdB, X_{1,peak}=%.3f', ...
      x_opt(1), x_opt(2), x_opt(3), f_opt(1), f_opt(2)));
yline(ax, 0, 'k--', '0 dB');
yline(ax, -40, 'r--', '-40 dB');
plot(ax, Om(valid), TF_dB(valid), 'b-', 'LineWidth', 1.5);
xlim(ax, [0.1, Omega_Start]);

%% =========================================================================
%  多目标代价函数
% =========================================================================
function f = multiobj_func(x, sysP_template, Omega_Start, Omega_Next, Nsteps)
    % 极大惩罚值 (代表求解失败或曲线发散)
    penalty = [1e6, 1e6, 1e6, 1e6];
    
    sigma = x(1); kap_e = x(2); kap_c = x(3);
    sysP = sysP_template;
    sysP(8) = kap_e; sysP(9) = kap_c; sysP(10) = sigma;
    global Fw
    
    try
        y_init = [zeros(15,1); Omega_Start];
        [x0_full, ok0] = newton('nondim_temp2', y_init, sysP);
        if ~ok0, f = penalty; return; end
        x0 = x0_full(1:15);
        
        y_init2 = [x0; Omega_Next];
        [x1_full, ok1] = newton('nondim_temp2', y_init2, sysP);
        if ~ok1, f = penalty; return; end
        x1 = x1_full(1:15);
        
        [x_res, ~] = branch_follow2('nondim_temp2', Nsteps, Omega_Start, Omega_Next, x0, x1, sysP);
    catch
        f = penalty; return;
    end
    
    Om = x_res(16,:).';
    
    % --- 计算传递力 TF_dB ---
    be2 = sysP(2); mu = sysP(3); ze2 = sysP(6); ga2 = sysP(11);
    x2 = x_res(6:10,:).';
    x2_dot = zeros(size(x2));
    x2_dot(:,2) = Om .* x2(:,3); x2_dot(:,3) = -Om .* x2(:,2);
    x2_dot(:,4) = 3*Om .* x2(:,5); x2_dot(:,5) = -3*Om .* x2(:,4);
    x2_cub = cubic_proj_013_batch(x2);
    ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
    ft_amp = hypot(hypot(ft(:,2), ft(:,3)), hypot(ft(:,4), ft(:,5)));
    TF_dB = 20*log10(max(ft_amp ./ Fw, 1e-300));
    
    % --- 计算上层位移 X1 ---
    x1_coeff = x_res(1:5,:).';
    X1_amp = hypot(hypot(x1_coeff(:,2), x1_coeff(:,3)), hypot(x1_coeff(:,4), x1_coeff(:,5)));
    
    valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
    Om_v = Om(valid); TF_v = TF_dB(valid); X1_v = X1_amp(valid);
    
    if isempty(Om_v)
        f = penalty; return;
    end
    
    % --- 目标 1: 峰值传递率 ---
    obj1 = max(TF_v);
    
    % --- 目标 2: 峰值动态位移 ---
    obj2 = max(X1_v);
    
    % --- 目标 3: 单位隔离频带起始频率 (越小代表高频隔振区越宽) ---
    idx_above_0 = find(TF_v >= 0);
    if isempty(idx_above_0)
        obj3 = min(Om_v); 
    else
        obj3 = max(Om_v(idx_above_0));
    end
    
    % --- 目标 4: -40dB 隔离频带起始频率 ---
    idx_above_m40 = find(TF_v >= -40);
    if isempty(idx_above_m40)
        obj4 = min(Om_v);
    else
        obj4 = max(Om_v(idx_above_m40));
    end
    
    f = [obj1, obj2, obj3, obj4];
end

%% ============ 辅助函数 ============
function cubic = cubic_proj_013_batch(U)
    [~, T_mat, T_inv] = get_AFT_matrices_local();
    cubic = (T_inv * ( (T_mat * U.').^3 )).';
end

function [N, T_mat, T_inv] = get_AFT_matrices_local()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64; t = (0:pN-1)'*(2*pi/pN);
        c1=cos(t); s1=sin(t); c3=cos(3*t); s3=sin(3*t); dc=ones(pN,1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv; pTinv(1,:) = (1/pN) * dc';
    end
    N = pN; T_mat = pT; T_inv = pTinv;
end