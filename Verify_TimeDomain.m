%% Verify_TimeDomain_Jump.m
% 时域纯数值积分验证：探究系统在悬崖边缘 (Omega=0.26, Fw=0.008) 的真实物理状态
% 包含：稳态时间历程图 (Time history) + 相轨迹与庞加莱截面 (Phase portrait & Poincare section)

clc; clear; close all;

%% -------- 1. 系统参数定义 (与参数严格对齐) --------
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.2;       % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 反推 v 与非线性系数
v = 2.5;        
alpha1 = v - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

% 模型系数
be1 = 1.0; 
al1 = alpha1 - be1; 
be2 = alpha2;
ga1 = gamma1; 
ga2 = gamma2; 
ze2 = 0.05;     % 下层阻尼比

% 电路参数 (带控制)
lam   = 0.18;
kap_e = 0.395;
kap_c = 0.032;
sigma = 0.623;

theta = sqrt(lam);

% 外部激励 (锁定在悬崖边缘的失稳测试点)
Omega = 2.64;
Fw    = 0.005;

%% -------- 2. 数值积分设置 --------
T_period = 2*pi/Omega;      % 激励周期
N_periods = 1000;           % 总仿真周期数（悬崖边缘瞬态过程长，必须给足时间）
t_span = [0, N_periods * T_period];

% 初始条件 [x1, v1, x2, v2, q, i_q] - 赋予微小初始扰动
y0 = [0.01; 0; 0; 0; 0; 0]; 

% 提高 ODE45 求解精度，防止非线性奇异导致数值耗散
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

fprintf('正在进行 ode45 数值积分 (Omega=%.2f, Fw=%.3f)...\n', Omega, Fw);
fprintf('总仿真时间: %.1f 个周期，请稍候...\n', N_periods);

% 求解系统微分方程
sol = ode45(@(t, y) sys_dynamics(t, y, be1, al1, be2, ga1, ga2, mu, ze2, theta, kap_e, kap_c, sigma, Omega, Fw), t_span, y0, options);

%% -------- 3. 提取稳态数据与计算庞加莱截面 --------
% 取最后 100 个周期作为绝对稳态数据
N_steady = 100;
t_steady_start = (N_periods - N_steady) * T_period;
t_steady_end   = N_periods * T_period;

% 构造密集的时间点用于画平滑的时间历程和相轨迹
t_dense = linspace(t_steady_start, t_steady_end, 10000);
y_dense = deval(sol, t_dense);

% 构造庞加莱截面 (Poincare Map) 的采样时间点：严格每隔一个周期 T 采样一次
t_poincare = t_steady_start : T_period : t_steady_end;
y_poincare = deval(sol, t_poincare);

%% -------- 4. SCI 风格联合绘图 --------
fontName = 'Times New Roman';
fsLab = 13; fsTit = 13; fsLeg = 11;

figure('Color', 'w', 'Position', [100, 100, 1000, 450]);
tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'loose');

% === 子图 1：稳态时间历程 (Time History) ===
nexttile; hold on; box on; grid on;
% 绘制 m1 和 m2 的位移
plot(t_dense / T_period, y_dense(1,:), 'b-', 'LineWidth', 1.2, 'DisplayName', 'Upper mass $x_1$');
plot(t_dense / T_period, y_dense(3,:), 'r--', 'LineWidth', 1.2, 'DisplayName', 'Lower mass $x_2$');
xlabel('Time / Excitation Periods ($t/T$)', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsLab);
ylabel('Displacement', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsLab);
title(sprintf('Steady-state Time History ($\\Omega=%.2f, F_w=%.3f$)', Omega, Fw), ...
    'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsTit);
legend('Location', 'best', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsLeg);
set(gca, 'FontName', fontName, 'FontSize', 11);
xlim([N_periods - 20, N_periods]); % 只展示最后20个周期以便看清波形

% === 子图 2：相轨迹与庞加莱截面 (Phase Portrait & Poincare) ===
nexttile; hold on; box on; grid on;
% 绘制 m1 的相轨迹 (连续实线)
plot(y_dense(1,:), y_dense(2,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8, 'DisplayName', 'Phase Trajectory');
% 叠加庞加莱截面点 (散点)，已移除转义字符防止报错
scatter(y_poincare(1,:), y_poincare(2,:), 25, 'ro', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Poincare Map');

xlabel('Displacement $x_1$', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsLab);
ylabel('Velocity $\dot{x}_1$', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsLab);
title('Phase Portrait & Poincare Section', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsTit);
legend('Location', 'best', 'Interpreter', 'latex', 'FontName', fontName, 'FontSize', fsLeg);
set(gca, 'FontName', fontName, 'FontSize', 11);

fprintf('计算与绘图完成！\n');

%% -------- 局部函数：状态空间微分方程 --------
function dy = sys_dynamics(t, y, be1, al1, be2, ga1, ga2, mu, ze2, theta, kap_e, kap_c, sigma, Omega, Fw)
    % 状态向量分配: y = [x1; v1; x2; v2; q; q_dot]
    x1 = y(1); v1 = y(2);
    x2 = y(3); v2 = y(4);
    q  = y(5); iq = y(6);
    
    dx = x1 - x2;
    dv = v1 - v2;
    
    % 机械恢复力
    f12 = (be1 + al1)*dx + ga1*(dx^3);
    f2g = be2*x2 + ga2*(x2^3) + 2*mu*ze2*v2;
    
    % 电磁力与电学方程
    dy = zeros(6,1);
    dy(1) = v1;
    dy(2) = -f12 + theta*iq + Fw*cos(Omega*t);                  % 上层加速度
    dy(3) = v2;
    dy(4) = (f12 - f2g - theta*iq) / mu;                        % 下层加速度
    dy(5) = iq;
    dy(6) = (-sigma*iq - kap_c*q - theta*dv) / kap_e;           % 电路响应 (保证 kap_e > 0)
end