%% 纯非线性测试
% 已移除失效的逐点扫频与 newton_enhanced
% 依赖外部函数: nondim_temp2, newton, branch_follow2

clc; clear; close all;

%% -------- 1. 基础参数与电路参数定义 --------
% Wang 图注给定参数（BG Model）
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0;     % 下层水平弹簧刚度比
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比

% 反推 v 与非线性系数
v = 2.5;        % 由 L=4/9, K1=1 反推
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
P.ze1 = 0.05;   % 下层阻尼比

% 待验证的电路参数（原第二组参数）
P.lam   = 0.0;
P.kap_e = 0.0;
P.kap_c = 0.0;
P.sigma = 0.0;

% 组装系统参数向量 sysP
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw
Fw = 0.005;

%% -------- 2. 弧长延拓起步设置 --------
Omega_Start = 10.0;
Omega_Step  = -0.01;  % 弧长法自适应步长，初始给定合理负增量即可
Omega_Next  = Omega_Start + Omega_Step;

fprintf('开始计算带电路参数单组曲线...\n');
fprintf('参数: lam=%.2f, kap_e=%.2f, kap_c=%.2f, sigma=%.2f\n', ...
        P.lam, P.kap_e, P.kap_c, P.sigma);

%% -------- 3. 求解初始两个点 --------
% 第一个点 (高频起点)
y_init = zeros(15,1);
y_init(end+1) = Omega_Start;
[x0_full, ok0] = newton('nondim_temp2', y_init, sysP);

if ~ok0
    error('高频起点求解失败，请检查参数或初始猜测值。');
end
x0 = x0_full(1:15);

% 第二个点 (用于确定弧长初始切线方向)
y_init2 = [x0; Omega_Next];
[x1_full, ok1] = newton('nondim_temp2', y_init2, sysP);

if ~ok1
    error('第二个初始点求解失败，无法启动弧长法。');
end
x1 = x1_full(1:15);

%% -------- 4. 弧长延拓主循环 --------
% 设定延拓步数为 3000 步，足以覆盖到极低频
[x_res, ~] = branch_follow2('nondim_temp2', 3000, Omega_Start, Omega_Next, x0, x1, sysP);

%% -------- 5. 计算力传递率 TF --------
Om  = x_res(16,:).';
be2 = sysP(2);
mu  = sysP(3);
ze2 = sysP(6);
ga2 = sysP(11);

% 提取下层位移的谐波系数 [x20, a21, b21, a23, b23]
x2 = x_res(6:10,:).';    

% 计算 x2_dot
W = Om;
x2_dot = zeros(size(x2));
x2_dot(:,1) = 0;
x2_dot(:,2) = W .* x2(:,3);
x2_dot(:,3) = -W .* x2(:,2);
x2_dot(:,4) = 3*W .* x2(:,5);
x2_dot(:,5) = -3*W .* x2(:,4);

% 批量计算三次非线性项投影
x2_cub = cubic_proj_013_batch(x2);

% 合成传递力 (基波与三次谐波)
ft = be2*x2 + ga2*x2_cub + 2*mu*ze2*x2_dot;
ft1 = hypot(ft(:,2), ft(:,3));
ft3 = hypot(ft(:,4), ft(:,5));
ft_amp = hypot(ft1, ft3);

% 转换为分贝
TF    = ft_amp ./ Fw;
TF_dB = 20*log10(max(TF, 1e-300));

% 剔除无效点（NaN、Inf 以及超出左边界的负频率点）
valid = isfinite(Om) & isfinite(TF_dB) & (Om > 0);
Om_valid = Om(valid);
TF_dB_valid = TF_dB(valid);

%% -------- 6. 绘图 --------
figure('Color','w', 'Position',[150 150 700 500]);
ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'XScale','log');
xlabel(ax, '\Omega (log scale)');
ylabel(ax, 'Force Transmissibility 20log_{10}(|f_t|/f) (dB)');
title_str = sprintf('BG Model (Electromechanical): \\lambda=%.2f, \\kappa_e=%.2f, \\kappa_c=%.2f, \\sigma=%.2f', ...
                    P.lam, P.kap_e, P.kap_c, P.sigma);
title(ax, title_str);

yline(ax, 0, 'k--', '0 dB');
% 绘制完整的频响曲线（包含可能存在的多值折叠区域）
plot(ax, Om_valid, TF_dB_valid, 'b-', 'LineWidth', 1.5);

% 设置视角边界限制以便于观察
xlim(ax, [0.1, Omega_Start]);
hold(ax,'off');

%% ============ 辅助函数：AFT 批量计算立方项 ============
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