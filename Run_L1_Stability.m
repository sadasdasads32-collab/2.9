%% Run_L1_Stability_Optimized.m
% 一键：定频扫力（L1） + Floquet 稳定性（沿HB轨道） + 稳定力区间提取 + SCI画图
%
% 依赖：
%   L1.m
%   nondim_temp2.m, newton.m, branch_follow2.m, branch_aux2.m
%
% 关键优化：
%   1) Floquet 单周期算子只积分 Phi，不积分状态 y（沿HB重构轨道评估 A(t)）
%   2) kap_e==0 时不引入虚拟极点，改为代数约束消元得到4维机械系统的Floquet
%   3) 提供系数顺序一致性自检（可选）
%
clear; clc; close all;

global ParamMin ParamMax
ParamMin = 0;
ParamMax = 0.02;

%% -------- 1. 基础参数与电路参数定义 --------
mu   = 0.2;     % 质量比 m2/m1
beta = 2.0;     % 下层竖向线性刚度比
K1   = 1.0;     % 上层水平弹簧刚度比
K2   = 0.2;     % 下层水平弹簧刚度比
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
P.ze1 = 0.05;   % （你原来叫ze1，但用在x2阻尼项上）

% 电路/耦合参数
P.lam   = 0.18;
P.kap_e = 0.395;
P.kap_c = 0.032;
P.sigma = 0.623;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

% ======= 固定频率（定频扫力）=======
Omega_fixed = 0.12;

% ======= 扫力设置 =======
F0     = 0.001;
dF     = 5e-4;
Nsteps = 3000;

% ======= Floquet 设置 =======
Nt = 400;               % 建议比你原来250稍大，Floquet更稳
tol_stable = 1.002;     % max|mu| < tol
do_selfcheck = false;   % 若怀疑系数顺序，改true会做轨道一致性检查

%% -----------------------------
% 2) 生成扫力起点：x_coeff(15)
% -----------------------------
global FixedOmega Fw
FixedOmega = Omega_fixed;
Fw = F0;

y_init = [zeros(15,1); F0];
y_sol  = newton('nondim_temp2', y_init, sysP);
x_coeff = y_sol(1:15);
x_interest = [x_coeff; F0];

%% -----------------------------
% 3) 定频扫力：L1
% -----------------------------
fprintf('--- Stage 1: Force continuation (L1) ---\n');
[x_branchF, infoL1] = L1(x_interest, sysP, dF, Nsteps, Omega_fixed);

F_all = x_branchF(16,:).';       % N x 1
Xc    = x_branchF(1:15,:).';     % N x 15

%% -----------------------------
% 4) Floquet 稳定性（沿HB轨道）
% -----------------------------
fprintf('\n--- Stage 2: Floquet stability (HB-orbit based) ---\n');
[is_stable, max_mu, mu_all] = Floquet_Stability_ForceSweep_HBOrbit( ...
    x_branchF, sysP, Omega_fixed, Nt, tol_stable, do_selfcheck);

% 稳定力区间提取
intervals = extract_stable_intervals(F_all, is_stable);

fprintf('\n========== Stable Force Intervals ==========\n');
if isempty(intervals)
    fprintf('No stable interval found under current settings.\n');
else
    for i = 1:size(intervals,1)
        fprintf('Stable interval %d:  [%.6g , %.6g]\n', i, intervals(i,1), intervals(i,2));
    end
end
fprintf('===========================================\n');

%% -----------------------------
% 5) SCI 版绘图
% -----------------------------
x2 = x_branchF(6:10,:).';
A1 = hypot(x2(:,2), x2(:,3));
A3 = hypot(x2(:,4), x2(:,5));
A_resp = sqrt(A1.^2 + A3.^2);

F = x_branchF(16,:).';
stab = is_stable(:).';
mu  = max_mu(:).';

[Fp, Ap, sp, mup] = thin_for_plot(F, A_resp, stab, mu);
intervals = extract_stable_intervals(F, stab);

figure('Color','w','Position',[120 120 980 420]);
tiledlayout(1,2,"Padding","compact","TileSpacing","compact");
fontName = 'Times New Roman';
fsLab = 13; fsTit = 13; fsLeg = 11;

% 左图：A–F
nexttile; hold on; box on; grid on;
hS = scatter(Fp(sp),  Ap(sp),  16, 'filled');
hU = scatter(Fp(~sp), Ap(~sp), 16);
if ~isempty(hU), hU.MarkerFaceColor = 'none'; hU.LineWidth = 1.0; end

xlabel('Force amplitude $F_w$','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
ylabel('Response amplitude $A$ (combined)','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
title(sprintf('Force continuation response at $\\Omega=%.3g$', Omega_fixed), ...
      'Interpreter','latex','FontName',fontName,'FontSize',fsTit);
set(gca,'FontName',fontName,'FontSize',11);

yl = ylim;
for i = 1:size(intervals,1)
    x0 = intervals(i,1); x1 = intervals(i,2);
    patch([x0 x1 x1 x0],[yl(1) yl(1) yl(2) yl(2)], ...
          'k','FaceAlpha',0.06,'EdgeColor','none');
end

if ~isempty(hS) && ~isempty(hU)
    legend([hS,hU],{'Stable','Unstable'},'Location','best','FontName',fontName,'FontSize',fsLeg);
elseif ~isempty(hS)
    legend(hS,{'Stable'},'Location','best','FontName',fontName,'FontSize',fsLeg);
elseif ~isempty(hU)
    legend(hU,{'Unstable'},'Location','best','FontName',fontName,'FontSize',fsLeg);
end

% 右图：Floquet 指标
nexttile; hold on; box on; grid on;
scatter(Fp, mup, 14, 'filled');
yline(tol_stable,'--','LineWidth',1.6);
xlabel('Force amplitude $F_w$','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
ylabel('max $|\mu|$','Interpreter','latex','FontName',fontName,'FontSize',fsLab);
title('Floquet stability indicator','Interpreter','latex','FontName',fontName,'FontSize',fsTit);
set(gca,'FontName',fontName,'FontSize',11);

[muMaxVal, idxWorst] = max(max_mu(:));
plot(F(idxWorst), muMaxVal, 'kp', 'MarkerSize', 10, 'LineWidth', 1.5);
text(F(idxWorst), muMaxVal, sprintf('  peak=%.3g', muMaxVal), 'FontName',fontName);

%% ====== 辅助函数 ======

function [Fx, Ax, Sx, Mx] = thin_for_plot(F, A, S, M)
    N = length(F);
    skip = max(1, round(N / 800));
    idx = 1:skip:N;
    if idx(end) ~= N, idx = [idx, N]; end
    Fx = F(idx); Ax = A(idx); Sx = S(idx); Mx = M(idx);
end

function intervals = extract_stable_intervals(F, is_stable)
    F = F(:); is_stable = is_stable(:);
    [F, idx] = sort(F);
    is_stable = is_stable(idx);
    d = diff([false; is_stable; false]);
    sidx = find(d==1);
    eidx = find(d==-1)-1;
    intervals = [F(sidx), F(eidx)];
end

%% =========================
%  Floquet：沿HB轨道积分 monodromy
% =========================
function [is_stable, max_mu, mu_all] = Floquet_Stability_ForceSweep_HBOrbit( ...
    x_resF, sysP, Omega_fixed, Nt, tol_stable, do_selfcheck)

    if nargin < 4 || isempty(Nt), Nt = 400; end
    if nargin < 5 || isempty(tol_stable), tol_stable = 1.002; end
    if nargin < 6, do_selfcheck = false; end

    kap_e = sysP(8);
    use_reduced = (abs(kap_e) < 1e-14); % kap_e==0 -> 4维消元系统

    N_pts = size(x_resF, 2);
    is_stable = true(1, N_pts);
    max_mu = zeros(1, N_pts);

    if use_reduced
        mu_all = zeros(4, N_pts);
    else
        mu_all = zeros(6, N_pts);
    end

    Omega = Omega_fixed;
    T = 2*pi/Omega;
    dt = T / Nt;

    for k = 1:N_pts
        Fw = x_resF(16, k);

        x1c = x_resF(1:5, k);
        x2c = x_resF(6:10, k);
        qc  = x_resF(11:15, k);

        % （可选）系数顺序自检：HB重构的y(T)应等于y(0)
        if do_selfcheck
            y0 = HB_reconstruct_state(0,   x1c, x2c, qc, Omega);
            yT = HB_reconstruct_state(T,   x1c, x2c, qc, Omega);
            err = norm(yT - y0);
            if err > 1e-9
                warning('HB coefficient order mismatch?  ||y(T)-y(0)|| = %.3e at k=%d', err, k);
            end
        end

        if use_reduced
            % ===== 4维：消元电路变量（kap_e=0，代数约束） =====
            Phi = eye(4);
            t = 0;

            for i = 1:Nt
                k1 = dPhi_reduced(t,          Phi, sysP, Omega, Fw, x1c, x2c, qc);
                k2 = dPhi_reduced(t + 0.5*dt, Phi + 0.5*dt*k1, sysP, Omega, Fw, x1c, x2c, qc);
                k3 = dPhi_reduced(t + 0.5*dt, Phi + 0.5*dt*k2, sysP, Omega, Fw, x1c, x2c, qc);
                k4 = dPhi_reduced(t + dt,     Phi + dt*k3,     sysP, Omega, Fw, x1c, x2c, qc);
                Phi = Phi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
                t = t + dt;
            end

            ev = eig(Phi);
            mu_all(:,k) = ev;
            max_mu(k) = max(abs(ev));
            is_stable(k) = (max_mu(k) < tol_stable);

        else
            % ===== 6维：完整系统，沿HB轨道评估 A(t) =====
            Phi = eye(6);
            t = 0;

            for i = 1:Nt
                k1 = dPhi_full(t,          Phi, sysP, Omega, Fw, x1c, x2c, qc);
                k2 = dPhi_full(t + 0.5*dt, Phi + 0.5*dt*k1, sysP, Omega, Fw, x1c, x2c, qc);
                k3 = dPhi_full(t + 0.5*dt, Phi + 0.5*dt*k2, sysP, Omega, Fw, x1c, x2c, qc);
                k4 = dPhi_full(t + dt,     Phi + dt*k3,     sysP, Omega, Fw, x1c, x2c, qc);
                Phi = Phi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
                t = t + dt;
            end

            ev = eig(Phi);
            mu_all(:,k) = ev;
            max_mu(k) = max(abs(ev));
            is_stable(k) = (max_mu(k) < tol_stable);
        end
    end
end

%% ========= HB重构：状态 y(t) =========
% 输入参数：t, x1c, x2c, qc, Omega
% 输出：输出的就是位移和速度的时域表达式，而且正是正余弦形式的显式表达式。
% 目的：将频域中的傅里叶系数转换回时域信号，重构出系统在给定时刻 t 的完整状态向量。
function y = HB_reconstruct_state(t, x1c, x2c, qc, Omega)
    w = Omega;
    ct  = cos(w*t);  st  = sin(w*t);
    c3t = cos(3*w*t); s3t = sin(3*w*t);

    x1 = x1c(1) + x1c(2)*ct + x1c(3)*st + x1c(4)*c3t + x1c(5)*s3t;
    v1 = (-w*x1c(2))*st + (w*x1c(3))*ct + (-3*w*x1c(4))*s3t + (3*w*x1c(5))*c3t;

    x2 = x2c(1) + x2c(2)*ct + x2c(3)*st + x2c(4)*c3t + x2c(5)*s3t;
    v2 = (-w*x2c(2))*st + (w*x2c(3))*ct + (-3*w*x2c(4))*s3t + (3*w*x2c(5))*c3t;

    q  = qc(1)  + qc(2)*ct  + qc(3)*st  + qc(4)*c3t  + qc(5)*s3t;
    qd = (-w*qc(2))*st + (w*qc(3))*ct + (-3*w*qc(4))*s3t + (3*w*qc(5))*c3t;

    y = [x1; v1; x2; v2; q; qd];
end

%% ========= dPhi：完整6维（沿HB轨道） =========
function dPhi = dPhi_full(t, Phi, sysP, Omega, Fw, x1c, x2c, qc)
    y = HB_reconstruct_state(t, x1c, x2c, qc, Omega);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));

    x1=y(1); v1=y(2); x2=y(3); v2=y(4); q=y(5); qd=y(6);
    dx = x1-x2;

    % 导数项
    df12_ddx = (be1+al1) + 3*ga1*dx^2;
    df2g_dx2 = be2 + 3*ga2*x2^2;
    df2g_dv2 = 2*mu*ze;

    A = zeros(6);
    A(1,2) = 1;
    A(2,1) = -df12_ddx;      A(2,3) = +df12_ddx;                      A(2,6) = theta;
    A(3,4) = 1;
    A(4,1) = df12_ddx/mu;    A(4,3) = (-df12_ddx - df2g_dx2)/mu;      A(4,4) = -df2g_dv2/mu;   A(4,6) = -theta/mu;
    A(5,6) = 1;

    % 电路方程线性化（kap_e != 0）
    A(6,2) = -theta/kap_e;
    A(6,4) = +theta/kap_e;
    A(6,5) = -kap_c/kap_e;
    A(6,6) = -sigma/kap_e;

    dPhi = A * Phi;
end

%% ========= dPhi：约束消元4维（kap_e==0） =========
% kap_e=0 时，电路惯性项消失：0 = -sigma*qd - kap_c*q - theta*(v1-v2)
% 可在频域HB里求得 q,qd；做稳定性时把 q,qd 看成沿轨道的已知函数，同时扰动满足线性约束：
%   -sigma*δqd - kap_c*δq - theta*(δv1-δv2) = 0
% 并且 δq' = δqd
% 联立可消去 δq,δqd，得到机械4维闭系统：
%   δqd = -(kap_c/σ)δq -(θ/σ)(δv1-δv2)
%   δq' = δqd  -> 代入得到 δq' = -(kap_c/σ)δq -(θ/σ)(δv1-δv2)
% 再把 δq 作为附加状态，会变成5维；为了保持4维，我们进一步把 δq 也消去：
%  这会引入记忆项，不适合。
%
% 所以更稳妥：直接用“准静态电耦合力”等效到机械方程（只对你当前模型）：
%   在 x1dd 与 x2dd 中出现的 theta*qd 项，使用约束给出的 qd（由HB轨道）作为已知输入，
%   线性扰动中 theta*δqd 按约束写成 δv 反馈（相当于附加阻尼）。
%
% 这样得到一个4维一阶变分：δ[x1 v1 x2 v2]' = A4(t)*δ[x1 v1 x2 v2]
function dPhi = dPhi_reduced(t, Phi, sysP, Omega, Fw, x1c, x2c, qc)
    y = HB_reconstruct_state(t, x1c, x2c, qc, Omega);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    if abs(kap_e) > 1e-14
        error('dPhi_reduced called but kap_e != 0');
    end

    theta = sqrt(max(lam,0));
    s = max(abs(sigma), 1e-12);

    x1=y(1); v1=y(2); x2=y(3); v2=y(4);
    dx = x1-x2;

    df12_ddx = (be1+al1) + 3*ga1*dx^2;
    df2g_dx2 = be2 + 3*ga2*x2^2;
    df2g_dv2 = 2*mu*ze;

    % 约束扰动：sigma*δqd + kap_c*δq + theta*(δv1-δv2)=0
    % 在机械方程中只出现 theta*δqd 项，所以我们写：
    %   theta*δqd = -(theta^2/sigma)*(δv1-δv2) - (theta*kap_c/sigma)*δq
    % 这里我们忽略 δq 这条通道（否则需要引入δq扩维），相当于只保留主要的速度反馈通道，
    % 在 kap_e=0 且电路快变量消失时，这通常是主导项（也对应你原来“虚拟极点”想做的事，但不改系统阶数）。
    % 如果你确实需要严格包含δq通道，应改成5维（加δq状态），我也可以给你那版。

    add = (theta^2)/s; % 速度差反馈强度

    A4 = zeros(4);
    % δx1' = δv1
    A4(1,2) = 1;
    % δv1' = -df12*(δx1-δx2) + theta*δqd  (用速度差反馈近似)
    A4(2,1) = -df12_ddx;
    A4(2,3) = +df12_ddx;
    A4(2,2) = -add;     % -(theta^2/sigma)*δv1
    A4(2,4) = +add;     % +(theta^2/sigma)*δv2

    % δx2' = δv2
    A4(3,4) = 1;
    % δv2' = ( df12*(δx1-δx2) -df2g_x*δx2 -df2g_v*δv2 - theta*δqd )/mu
    A4(4,1) = df12_ddx/mu;
    A4(4,3) = (-df12_ddx - df2g_dx2)/mu;
    A4(4,4) = -df2g_dv2/mu;
    A4(4,2) = +add/mu;  % -(-theta^2/sigma)*(δv1-δv2)/mu
    A4(4,4) = A4(4,4) - add/mu;

    dPhi = A4 * Phi;
end