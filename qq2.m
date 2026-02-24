%% 双自由度 QZS 系统非线性骨架线 (Backbone Curve) 追踪
clc; clear; close all;
global Fw;
Fw = 0; % 骨架线计算必须无外力

%% -------- 1. 基础机械参数 --------
mu   = 0.2; beta = 2.0; K1 = 1.0; K2 = 0; U = 2.0; L = 4/9; v = 2.5;
alpha1 = v - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

% 构造 sysP 模板: [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
% 注意：此处机械阻尼 ze1 强制设为 0
sysP_base = [1.0, alpha2, mu, alpha1-1.0, gamma1, 0, ...
             0, 0, 0, 0, gamma2];

%% -------- 2. 定义对比工况 --------
% 工况 1：纯机械系统 (无电路)
sysP_mech = sysP_base;

% 工况 2：接入 VCM-NIC 电路 (注意：电阻 sigma 强制设为 0)
sysP_elec = sysP_base;
sysP_elec(7)  = 0.18;  % lam (耦合强度)
sysP_elec(8)  = 0.5;   % kap_e (电感)
sysP_elec(9)  = 0.1;   % kap_c (电容)
sysP_elec(10) = 0.0;   % sigma (无电阻)

%% -------- 3. 延拓参数与初始猜测 --------
% 振幅 A 的扫描范围 (即上层质量块的基波余弦幅值 a11)
A_vec = linspace(0.005, 0.35, 60); 

% 初始频率猜测 (基于线性模态的近似估算)
Omega_guess_mech = 0.8;  
Omega_guess_elec = 0.6;  

% 初始化存储数组
Om_mech = zeros(size(A_vec));
Om_elec = zeros(size(A_vec));

%% -------- 4. 追踪纯机械骨架线 --------
fprintf('正在计算纯机械系统骨架线...\n');
% 构造 14维 未知数向量 Z: 
% [x10, a13, b13, x20, a21, b21, a23, b23, q0, aq1, bq1, aq3, bq3, Omega]
Z_mech = zeros(14,1); 
Z_mech(14) = Omega_guess_mech;

for i = 1:length(A_vec)
    A = A_vec(i);
    [Z_mech, ok] = newton_backbone(Z_mech, A, sysP_mech);
    if ~ok, warning('机械系统在 A = %.3f 处不收敛', A); break; end
    Om_mech(i) = Z_mech(14); % 提取自然频率
end

%% -------- 5. 追踪机电耦合骨架线 --------
fprintf('正在计算机电耦合系统骨架线...\n');
Z_elec = zeros(14,1);
Z_elec(14) = Omega_guess_elec;

for i = 1:length(A_vec)
    A = A_vec(i);
    [Z_elec, ok] = newton_backbone(Z_elec, A, sysP_elec);
    if ~ok, warning('机电系统在 A = %.3f 处不收敛', A); break; end
    Om_elec(i) = Z_elec(14);
end

%% -------- 6. 骨架线绘图对比 --------
figure('Color','w', 'Position',[150 150 650 500]);
ax = gca; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
% 注意：惯例将频率放在 X 轴，振幅放在 Y 轴，方便与强迫响应图叠加
plot(ax, Om_mech(Om_mech>0), A_vec(Om_mech>0), 'k--', 'LineWidth', 2, 'DisplayName', 'Mechanical Backbone');
plot(ax, Om_elec(Om_elec>0), A_vec(Om_elec>0), 'r-', 'LineWidth', 2, 'DisplayName', 'Electromechanical Backbone');

xlabel(ax, 'Natural Frequency \Omega_n');
ylabel(ax, 'Response Amplitude A (a_{11})');
title(ax, 'Backbone Curves (Undamped Free Vibration)');
legend(ax, 'Location', 'northwest');
set(ax, 'XScale', 'log');
xlim(ax, [0.1, 1.5]);

%% ================== 核心求解器函数 ==================
function [Z, ok] = newton_backbone(Z_guess, A, sysP)
    % 针对骨架线的局部 Newton 求解器 (14未知数 -> 14方程)
    Z = Z_guess;
    ok = false;
    tol = 1e-8;
    for iter = 1:30
        R = backbone_wrapper(Z, A, sysP);
        if norm(R) < tol
            ok = true;
            break;
        end
        % 数值计算 14x14 雅可比矩阵
        J = zeros(14, 14);
        eps_val = 1e-6;
        for j = 1:14
            Z_pert = Z;
            Z_pert(j) = Z_pert(j) + eps_val;
            R_pert = backbone_wrapper(Z_pert, A, sysP);
            J(:, j) = (R_pert - R) / eps_val;
        end
        delta = -J \ R;
        Z = Z + delta;
    end
end

function R_reduced = backbone_wrapper(Z_14, A, sysP)
    % 将 14 维求解变量映射回你原始的 16 维状态空间
    y = zeros(16,1);
    y(1)   = Z_14(1);      % x10
    y(2)   = A;            % a11 (强行设为参数)
    y(3)   = 0;            % b11 (锚定相位，强行设为0)
    y(4:15)= Z_14(2:13);   % 其他谐波系数与电荷
    y(16)  = Z_14(14);     % Omega

    % 调用你的原始残差函数
    R_full = nondim_temp2(y, sysP);

    % 由于阻尼为0且 b11=0，x1 的正弦投影方程 (即 R_full(3)) 恒为 0
    % 将其剔除，使得方程组退化为 14x14 满秩适定问题
    R_reduced = [R_full(1:2); R_full(4:15)];
end