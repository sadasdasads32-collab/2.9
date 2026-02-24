%% Compare_Circuit_Params.m
% 对比两组电路参数下的定频扫力响应及稳定性
% 参数组1: lam=0.18, kap_e=1.83, kap_c=0.01, sigma=0.54
% 参数组2: lam=0,   kap_e=0,   kap_c=0,   sigma=0
% 其余机械参数固定（同 Run_L1_Stability.m）
% 依赖：nondim_temp2.m, newton.m, L1.m, branch_follow2.m 等

clear; clc; close all;
global ParamMin ParamMax FixedOmega Fw
ParamMin = 0; ParamMax = 0.02;

%% -------- 固定机械参数（Wang 图15曲线1）--------
mu   = 0.2;     % 质量比
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0;       % 下层水平弹簧刚度比（原脚本 K2=0）
U    = 2.0;     % 几何非线性尺度参数
L    = 4/9;     % QZS 长度比
v    = 2.5;     % 调谐参数（由 L=4/9 反推）

alpha1 = v - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

P_base.be1 = 1.0;
P_base.al1 = alpha1 - P_base.be1;
P_base.be2 = alpha2;
P_base.ga1 = gamma1;
P_base.ga2 = gamma2;
P_base.mu  = mu;
P_base.ze1 = 0.05;   % 下层阻尼比

%% -------- 两组电路参数 --------
circuit_params = [
    0.18, 1.83, 0.01, 0.54;   % 组1：有电路
    0.0,  0.0,  0.0,  0.0     % 组2：无电路
    ];
nGroups = size(circuit_params, 1);
groupNames = {'With circuit (lam=0.18, kap_e=1.83, kap_c=0.01, sigma=0.54)', ...
              'No circuit (all zero)'};

%% -------- 扫力设置（固定频率）--------
Omega_fixed = 0.28;      % 固定频率（建议选在峰值附近）
F0     = 0.001;          % 初始力幅
dF     = 5e-4;           % 扫力步长（正值为增加，负值为减小；此处用正，分支会向右延伸）
Nsteps = 3000;           % 最大延拓步数

% Floquet 稳定性设置
Nt = 250;                % 每周期 RK4 步数
tol_stable = 1.002;      % 稳定判据：最大 Floquet 乘子模 < tol_stable

%% -------- 预存储结果 --------
results = cell(nGroups, 1);

%% -------- 循环处理两组参数 --------
for g = 1:nGroups
    fprintf('\n========== Processing group %d: %s ==========\n', g, groupNames{g});
    
    % 组装 sysP 向量（顺序必须与 nondim_temp2 一致）
    P = P_base;
    P.lam   = circuit_params(g,1);
    P.kap_e = circuit_params(g,2);
    P.kap_c = circuit_params(g,3);
    P.sigma = circuit_params(g,4);
    
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
            P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];
    
    % 生成扫力起点（在 F0 处求一个解）
    FixedOmega = Omega_fixed;
    Fw = F0;
    y_init = [zeros(15,1); F0];
    y_sol = newton('nondim_temp2', y_init, sysP);
    x_coeff = y_sol(1:15);
    x_interest = [x_coeff; F0];
    
    % 定频扫力延拓（调用 L1 函数）
    [x_branchF, infoL1] = L1(x_interest, sysP, dF, Nsteps, Omega_fixed);
    F_all = x_branchF(16,:).';
    Xc = x_branchF(1:15,:).';
    
    % Floquet 稳定性计算
    [is_stable, max_mu, mu_all] = Get_Floquet_Stability_ForceSweep( ...
        x_branchF, sysP, Omega_fixed, Nt, tol_stable);
    
    % 提取稳定力区间
    intervals = extract_stable_intervals(F_all, is_stable);
    
    % 计算响应幅值（以 x2 的基波+三次谐波合成幅值为例）
    x2 = x_branchF(6:10,:).';
    A1 = hypot(x2(:,2), x2(:,3));
    A3 = hypot(x2(:,4), x2(:,5));
    A_resp = sqrt(A1.^2 + A3.^2);
    
    % 存储结果
    results{g}.F = F_all;
    results{g}.A = A_resp;
    results{g}.is_stable = is_stable(:);
    results{g}.max_mu = max_mu(:);
    results{g}.intervals = intervals;
    
    % 显示稳定区间
    fprintf('Stable intervals for %s:\n', groupNames{g});
    if isempty(intervals)
        fprintf('  None\n');
    else
        for i = 1:size(intervals,1)
            fprintf('  [%.6g, %.6g]\n', intervals(i,1), intervals(i,2));
        end
    end
end

%% -------- 绘图对比 --------
% 准备绘图数据（抽稀，避免过密）
% 统一 F 范围
F_min = min(cellfun(@(r) min(r.F), results));
F_max = max(cellfun(@(r) max(r.F), results));

figure('Color','w','Position',[120 120 1200 500]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

fontName = 'Times New Roman';
fsLab = 13; fsTit = 13; fsLeg = 11;
colors = lines(2);
markers = {'o','s'};   % 两种标记，方便黑白打印区分

% 左图：A–F 响应对比
nexttile; hold on; box on; grid on;
for g = 1:nGroups
    F = results{g}.F;
    A = results{g}.A;
    stab = results{g}.is_stable;
    % 抽稀（仅用于绘图）
    [Fp, Ap, sp] = thin_for_plot(F, A, stab, []);
    
    % 稳定点：实心；不稳定点：空心
    scatter(Fp(sp), Ap(sp), 24, colors(g,:), 'filled', 'Marker', markers{g});
    scatter(Fp(~sp), Ap(~sp), 24, colors(g,:), 'Marker', markers{g}, ...
            'MarkerFaceColor', 'none', 'LineWidth', 1.0);
end
xlabel('Force amplitude $F_w$','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
ylabel('Response amplitude $A$ (combined)','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
title(sprintf('Force continuation response at $\\Omega=%.3g$', Omega_fixed), ...
      'Interpreter','latex','FontName',fontName,'FontSize',fsTit);
% 标注稳定区间（半透明色块）
yl = ylim;
for g = 1:nGroups
    intervals = results{g}.intervals;
    for i = 1:size(intervals,1)
        x0 = intervals(i,1); x1 = intervals(i,2);
        patch([x0 x1 x1 x0],[yl(1) yl(1) yl(2) yl(2)], ...
              colors(g,:), 'FaceAlpha',0.05, 'EdgeColor','none');
    end
end
uistack(findobj(gca,'Type','scatter'),'top'); % 散点置顶
% 图例（用稳定点样式代表整组）
h1 = scatter(nan, nan, 24, colors(1,:), 'filled', 'Marker', markers{1}, 'DisplayName', groupNames{1});
h2 = scatter(nan, nan, 24, colors(2,:), 'filled', 'Marker', markers{2}, 'DisplayName', groupNames{2});
legend([h1, h2], 'Location','best', 'FontName',fontName, 'FontSize',fsLeg);

% 右图：Floquet 指标对比
nexttile; hold on; box on; grid on;
for g = 1:nGroups
    F = results{g}.F;
    mu = results{g}.max_mu;
    [Fp, mup] = thin_for_plot(F, mu, [], []); % 抽稀
    plot(Fp, mup, '.-', 'Color', colors(g,:), 'Marker', markers{g}, 'LineWidth', 1.2, 'MarkerSize', 8);
end
yline(tol_stable, 'k--', 'LineWidth', 1.6);
xlabel('Force amplitude $F_w$','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
ylabel('max $|\mu|$','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
title('Floquet stability indicator','Interpreter','latex','FontName',fontName,'FontSize',fsTit);
set(gca,'FontName',fontName,'FontSize',11);
legend([h1, h2], 'Location','best', 'FontName',fontName, 'FontSize',fsLeg); % 复用左图图例句柄

% 统一 x 轴范围
xlim([F_min, F_max]);

%% ========== 辅助函数 ==========
function [Fx, Ax, Sx, Mx] = thin_for_plot(F, A, S, M)
    % 对数组进行间隔抽稀，最多保留约 800 点，同时保留端点
    N = length(F);
    skip = max(1, round(N / 800));
    idx = 1:skip:N;
    if idx(end) ~= N
        idx = [idx, N];
    end
    Fx = F(idx);
    Ax = A(idx);
    if ~isempty(S)
        Sx = S(idx);
    else
        Sx = [];
    end
    if ~isempty(M)
        Mx = M(idx);
    else
        Mx = [];
    end
end

function intervals = extract_stable_intervals(F, is_stable)
    F = F(:);
    is_stable = is_stable(:);
    [F, idx] = sort(F);
    is_stable = is_stable(idx);
    d = diff([false; is_stable; false]);
    sidx = find(d==1);
    eidx = find(d==-1)-1;
    intervals = [F(sidx), F(eidx)];
end

function [is_stable, max_mu, mu_all] = Get_Floquet_Stability_ForceSweep(x_resF, sysP, Omega_fixed, Nt, tol_stable)
    % 固定频率 Omega_fixed，对扫力分支 x_resF(16,k)=F 逐点计算 Floquet 乘子
    if nargin < 4 || isempty(Nt), Nt = 250; end
    if nargin < 5 || isempty(tol_stable), tol_stable = 1.002; end

    N_pts = size(x_resF, 2);
    is_stable = true(1, N_pts);
    max_mu = zeros(1, N_pts);
    mu_all = zeros(6, N_pts);

    for k = 1:N_pts
        Fw = x_resF(16, k);
        Omega = Omega_fixed;

        % 从HBM系数重构 t=0 初值
        x1c = x_resF(1:5, k); x2c = x_resF(6:10, k); qc = x_resF(11:15, k);

        x1_0 = x1c(1) + x1c(2) + x1c(4);
        v1_0 = Omega * x1c(3) + 3 * Omega * x1c(5);

        x2_0 = x2c(1) + x2c(2) + x2c(4);
        v2_0 = Omega * x2c(3) + 3 * Omega * x2c(5);

        q_0  = qc(1)  + qc(2)  + qc(4);
        qd_0 = Omega * qc(3)  + 3 * Omega * qc(5);

        y0 = [x1_0; v1_0; x2_0; v2_0; q_0; qd_0];

        % 变分矩阵初值
        Phi0 = eye(6);
        Y = [y0; Phi0(:)];

        % RK4 积分一周期
        T = 2*pi/Omega;
        dt = T/Nt;
        t = 0;
        for i = 1:Nt
            k1 = ext_ode(t,          Y,             sysP, Omega, Fw);
            k2 = ext_ode(t + 0.5*dt, Y + 0.5*dt*k1, sysP, Omega, Fw);
            k3 = ext_ode(t + 0.5*dt, Y + 0.5*dt*k2, sysP, Omega, Fw);
            k4 = ext_ode(t + dt,     Y + dt*k3,     sysP, Omega, Fw);
            Y = Y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
            t = t + dt;
        end

        PhiT = reshape(Y(7:end), 6, 6);
        ev = eig(PhiT);

        mu_all(:,k) = ev;
        max_mu(k) = max(abs(ev));
        is_stable(k) = (max_mu(k) < tol_stable);
    end
end

function dY = ext_ode(t, Y, sysP, Omega, Fw)
    % 扩展的 ODE：状态 + 变分矩阵
    y = Y(1:6);
    Phi = reshape(Y(7:end), 6, 6);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);
    theta = sqrt(max(lam,0));

    x1=y(1); v1=y(2); x2=y(3); v2=y(4); q=y(5); qd=y(6);
    dx = x1-x2; dv = v1-v2;

    f12 = (be1+al1)*dx + ga1*dx^3;
    f2g = be2*x2 + ga2*x2^3 + 2*mu*ze2*v2;

    x1dd = -f12 + theta*qd + Fw*cos(Omega*t);
    x2dd = ( f12 - f2g - theta*qd )/mu;

    if kap_e == 0
        tiny = 1e-12;
        qd_new = (-kap_c*q - theta*dv)/max(abs(sigma), tiny);
        qdd = (qd_new - qd)*50;
    else
        qdd = (-sigma*qd - kap_c*q - theta*dv)/kap_e;
    end

    dy = [v1; x1dd; v2; x2dd; qd; qdd];

    % Jacobian A(t)
    df12_ddx = (be1+al1) + 3*ga1*dx^2;
    df2g_dx2 = be2 + 3*ga2*x2^2;
    df2g_dv2 = 2*mu*ze2;

    A = zeros(6);
    A(1,2) = 1;
    A(2,1) = -df12_ddx;      A(2,3) = +df12_ddx;                      A(2,6) = theta;
    A(3,4) = 1;
    A(4,1) = df12_ddx/mu;    A(4,3) = (-df12_ddx - df2g_dx2)/mu;      A(4,4) = -df2g_dv2/mu;   A(4,6) = -theta/mu;
    A(5,6) = 1;

    if kap_e ~= 0
        A(6,2) = -theta/kap_e;   A(6,4) = theta/kap_e;    A(6,5) = -kap_c/kap_e;   A(6,6) = -sigma/kap_e;
    else
        s = max(abs(sigma),1e-12);
        A(6,2) = -50*theta/s;    A(6,4) = 50*theta/s;     A(6,5) = -50*kap_c/s;    A(6,6) = -50;
    end

    dPhi = A*Phi;
    dY = [dy; dPhi(:)];
end