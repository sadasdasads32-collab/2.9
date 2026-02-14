%% 非线性验证 (倒序扫频 + 力传递率)
% 目标：倒序扫频 (Omega: 10 -> 0.2)，绘制力传递率 TF = |Ft|/Fw 的 dB 曲线
% 横坐标：log scale (10 的幂次)

clear; clc; close all;

%% 1) 参数设置（你给的非线性测试参数：此处 ga1=0 等价线性）
P.be1 = 1;        % 上层基准刚度比（通常固定为 1）
P.be2 = 0.5;      % 下层对地等效刚度比
P.mu  = 0.2;      % 质量比 m2/m1

P.al1 = -0.5;     % 上层 QZS 的等效线性修正项（此处你仍保留 -0.9）
P.ga1 = 0.4883;        % 上层三次项（此处为 0：线性）
P.ga2 = 0.4883;      % 下层三次项（线性）

P.ze1 = 0.111803; % 下层阻尼比（你当前的定义体系）
P.lam   = 0.18;    
P.kap_e = 0.330426;
P.kap_c = 0.292043;
P.sigma = 0.250000;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.02;         % 外激励幅值（你给的）
FixedOmega = [];   % 扫频模式

%% 2) 倒序扫频范围
Omega_Start = 5.0;
Omega_End   = 0.2;
Omega_Step  = -0.01;   % 负步长

fprintf('开始倒序扫频: %.2f -> %.2f\n', Omega_Start, Omega_End);

%% 3) 高频起点求解（Newton）
y_init = zeros(15,1);
y_init(end+1) = Omega_Start;  % 16维：[15维HB状态 + Omega]

[x0_full, ok] = newton('nondim_temp2', y_init, sysP);
if ~ok, error('高频起点求解失败！'); end
x0 = x0_full(1:15);

%% 4) 构造第二个点（启动弧长法）
Omega_Next = Omega_Start + Omega_Step;
y_init2 = [x0; Omega_Next];

[x1_full, ok] = newton('nondim_temp2', y_init2, sysP);
if ~ok, error('第二个点求解失败！'); end
x1 = x1_full(1:15);

%% 5) 弧长延拓（倒序扫频）
[x_res, ~] = branch_follow2('nondim_temp2', 5000, Omega_Start, Omega_Next, x0, x1, sysP);

%% 6) 力传递率 TF（基于对地反力）+ dB + log横坐标
Om = x_res(16,:);

% HB 基波位移幅值 |X2|
% x2 在 6:10，基波 cos/sin 系数对应索引 7,8
A2 = sqrt(x_res(7,:).^2 + x_res(8,:).^2);

% 取参数
be2 = sysP(2);
mu  = sysP(3);
ze2 = sysP(6);   % 你的命名 ze1，但这里物理上对应 zeta2（下层对地阻尼比）

% 传递到地基的反力基波幅值 |Ft|
Ft_amp = A2 .* sqrt( be2^2 + (2*mu*ze2.*Om).^2 );

% 力传递率
TF = Ft_amp ./ Fw;
TF_dB = 20*log10(TF + 1e-300);

% 绘图：横坐标 10 的幂次（log scale）
figure('Color','w');
semilogx(Om, TF_dB, 'LineWidth', 1.8);
yline(0,'--'); grid on;

xlim([Omega_End Omega_Start]);
xlabel('\Omega (log scale)');
ylabel('Force Transmissibility 20log_{10}(|F_t|/F_w) (dB)');
title('Backward Sweep: Force Transmissibility');

% 可选：自定义刻度更像论文图（可删）
ax = gca;
ax.XTick = [0.1 0.2 0.5 1 2 5 10];
ax.XTickLabel = {'10^{-1}','2\times10^{-1}','5\times10^{-1}','10^{0}','2\times10^{0}','5\times10^{0}','10^{1}'};
