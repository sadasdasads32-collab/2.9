%% Run_L1_Floquet_vs_TimeDomain_ULTIMATE.m
% 目的：在同一条 A(F) 曲线上，把 Floquet 稳定/不稳定 与 时域跳跃区间对齐
% 终极版特性：
% 1) A(F) 统一使用 p2p 幅值（与准静态时域一致）
% 2) Floquet 指标图：散点 + 按F排序线（避免“回环柱子”视觉误导）
% 3) 对 L1 点做按F分箱去重（binning），抑制 fold 附近点堆叠造成的“橙色大柱子”
% 4) 自动提取 Floquet 不稳定窗口，并用淡灰阴影标注（两张图都标）
%
% 必需变量(workspace)：
%   x_branchF   : 16xN (L1输出，前15行HBM系数，最后一行为F)
%   Omega_fixed : 固定频率标量
%   sysP        : 11参数向量
% 可选变量(workspace)：
%   F_fwd,A_fwd,F_bwd,A_bwd : 时域准静态滞回曲线（建议叠加）
%
% 你只需要运行本脚本，不需要任何外部函数文件。



%% ========== 0) 检查输入 ==========
assert(exist('x_branchF','var')==1, 'Need x_branchF in workspace (from L1).');
assert(exist('Omega_fixed','var')==1, 'Need Omega_fixed in workspace.');
assert(exist('sysP','var')==1, 'Need sysP in workspace.');

F_L1 = x_branchF(16,:).';
N = numel(F_L1);
fprintf('L1 branch points: N=%d\n', N);

has_td = exist('F_fwd','var')==1 && exist('A_fwd','var')==1 && ...
         exist('F_bwd','var')==1 && exist('A_bwd','var')==1;

%% ========== 1) L1 幅值：用 p2p（与时域一致） ==========
x2c = x_branchF(6:10,:).';  % [a0 a1 b1 a3 b3]
A_L1_p2p = hb_p2p_from_x2coeff(x2c, Omega_fixed, 2048);

%% ========== 2) Floquet：沿 L1 分支逐点计算 max|mu| ==========
Nt = 800;            % 推荐 500~1200
tol_stable = 1.002;  % 稳定判据阈值（考虑数值误差）

fprintf('Computing Floquet multipliers along L1 branch... Nt=%d\n', Nt);
[is_stable_L1, max_mu_L1] = Floquet_along_L1(x_branchF, sysP, Omega_fixed, Nt, tol_stable);
fprintf('Floquet done.\n');

%% ========== 3) 按 F 分箱去重：抑制 fold 附近点堆叠 ==========
% 这一步是你图里“橙色大柱子”的核心修复点
dF_bin = 2e-5;  % 分箱宽度：1e-5 ~ 1e-4 可调（越小越保留细节）
[Fb, Ab, mub, stb] = bin_unique_by_F(F_L1, A_L1_p2p, max_mu_L1, is_stable_L1, dF_bin);

F_plot  = Fb;
A_plot  = Ab;
mu_plot = mub;
st_plot = stb;

fprintf('After binning: N=%d (from %d)\n', numel(F_plot), N);

%% ========== 4) 时域滞回 jump 点估计（可选） ==========
if has_td
    [Fup_est, Fdown_est] = estimate_jump_from_hysteresis(F_fwd,A_fwd,F_bwd,A_bwd);
    fprintf('Estimated jump from TD: F_down≈%.6g, F_up≈%.6g\n', Fdown_est, Fup_est);
end

%% ========== 5) Floquet 不稳定窗口提取（按F排序后识别） ==========
unstable_intervals = get_intervals_from_stability(F_plot, st_plot, false);

if ~isempty(unstable_intervals)
    fprintf('\n========== Floquet-UNSTABLE Force Intervals (binned) ==========\n');
    for i=1:size(unstable_intervals,1)
        fprintf('Unstable interval %d: [%.6g , %.6g]\n', i, unstable_intervals(i,1), unstable_intervals(i,2));
    end
else
    fprintf('\nNo unstable interval detected under tol_stable=%.4g (binned data)\n', tol_stable);
end

%% ========== 6) SCI 绘图 ==========
figure('Color','w','Position',[90 70 1120 540]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% ---- (a) A(F) with Floquet + TD ----
ax1 = nexttile; hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');

% 先画不稳定阴影（底层）
shade_intervals(ax1, unstable_intervals);

% L1 稳定/不稳定分段（按延拓顺序分段；这里用去重后的点）
plot_stable_unstable(ax1, F_plot, A_plot, st_plot, 2.8);

% 叠加时域准静态滞回
if has_td
    plot(ax1, F_fwd, A_fwd, 'LineWidth', 2.2);
    plot(ax1, F_bwd, A_bwd, 'LineWidth', 2.2);

    xline(ax1, Fdown_est, '--', sprintf('F_{down}=%.4g',Fdown_est), 'LineWidth', 1.2);
    xline(ax1, Fup_est,   '--', sprintf('F_{up}=%.4g',Fup_est),     'LineWidth', 1.2);
end

xlabel(ax1,'Force amplitude $F_w$','Interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel(ax1,'Response amplitude $A$ (p2p)','Interpreter','latex','FontSize',13,'FontName','Times New Roman');
title(ax1, sprintf('$A(F)$ with Floquet stability at $\\Omega=%.3g$', Omega_fixed), ...
      'Interpreter','latex','FontSize',13);

% 固定图例顺序（占位符）
h1 = plot(ax1, NaN,NaN,'-','LineWidth',2.8);
h2 = plot(ax1, NaN,NaN,'--','LineWidth',2.8);
h3 = plot(ax1, NaN,NaN,'-','LineWidth',2.2);
h4 = plot(ax1, NaN,NaN,'-','LineWidth',2.2);
legend(ax1,[h1,h2,h3,h4],{'L1 stable','L1 unstable','TD forward','TD backward'}, ...
    'Location','best','FontSize',11,'FontName','Times New Roman');

% ---- (b) Floquet indicator ----
ax2 = nexttile; hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');

% 不稳定阴影（底层）
shade_intervals(ax2, unstable_intervals);

% 散点：真实点云
scatter(ax2, F_plot, mu_plot, 14, 'filled');

% 排序线：视觉友好，不会“回环柱子”
[Fs, idx] = sort(F_plot);
plot(ax2, Fs, mu_plot(idx), 'LineWidth', 2.2);

yline(ax2, 1.0, '--', 'LineWidth', 1.2);

if has_td
    xline(ax2, Fdown_est, '--', 'LineWidth', 1.1);
    xline(ax2, Fup_est,   '--', 'LineWidth', 1.1);
end

xlabel(ax2,'Force amplitude $F_w$','Interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel(ax2,'$\max|\mu|$','Interpreter','latex','FontSize',13,'FontName','Times New Roman');
title(ax2,'Floquet stability indicator','FontSize',13,'FontName','Times New Roman');

%% ===================== Local functions =====================

function A_p2p = hb_p2p_from_x2coeff(x2c, Omega, Ns)
% 从HBM系数重构 x2(t)，输出 p2p 幅值
% x2c: Nx5 [a0 a1 b1 a3 b3]
    if nargin<3, Ns=2048; end
    t = linspace(0, 2*pi/Omega, Ns).';
    ct1 = cos(Omega*t); st1 = sin(Omega*t);
    ct3 = cos(3*Omega*t); st3 = sin(3*Omega*t);

    A_p2p = zeros(size(x2c,1),1);
    for k=1:size(x2c,1)
        a0=x2c(k,1); a1=x2c(k,2); b1=x2c(k,3); a3=x2c(k,4); b3=x2c(k,5);
        x = a0 + a1*ct1 + b1*st1 + a3*ct3 + b3*st3;
        A_p2p(k) = max(x)-min(x);
    end
end

function plot_stable_unstable(ax, x, y, is_stable, lw)
% 分段画稳定(实线)/不稳定(虚线)
    idx_change = find(diff(is_stable)~=0);
    seg_start = 1;
    for i=1:length(idx_change)+1
        if i>length(idx_change), seg_end=numel(x);
        else, seg_end=idx_change(i);
        end
        if is_stable(seg_start)
            plot(ax, x(seg_start:seg_end), y(seg_start:seg_end), '-', 'LineWidth', lw);
        else
            plot(ax, x(seg_start:seg_end), y(seg_start:seg_end), '--', 'LineWidth', lw);
        end
        seg_start = seg_end;
    end
end

function [is_stable, max_mu] = Floquet_along_L1(x_branchF, sysP, Omega_fixed, Nt, tol_stable)
% 沿 L1 分支逐点计算 Floquet 最大乘子模
    Npts = size(x_branchF,2);
    is_stable = true(Npts,1);
    max_mu = zeros(Npts,1);

    for k=1:Npts
        Fw = x_branchF(16,k);

        x1c = x_branchF(1:5,k); x2c = x_branchF(6:10,k); qc = x_branchF(11:15,k);
        Om = Omega_fixed;

        x1_0 = x1c(1)+x1c(2)+x1c(4);
        v1_0 = Om*x1c(3) + 3*Om*x1c(5);
        x2_0 = x2c(1)+x2c(2)+x2c(4);
        v2_0 = Om*x2c(3) + 3*Om*x2c(5);
        q_0  = qc(1)+qc(2)+qc(4);
        qd_0 = Om*qc(3) + 3*Om*qc(5);

        y0 = [x1_0; v1_0; x2_0; v2_0; q_0; qd_0];

        Phi0 = eye(6);
        Y = [y0; Phi0(:)];

        T = 2*pi/Om;
        dt = T/Nt;
        t = 0;

        for i=1:Nt
            k1 = ext_ode(t,          Y,             sysP, Om, Fw);
            k2 = ext_ode(t+0.5*dt,   Y+0.5*dt*k1,    sysP, Om, Fw);
            k3 = ext_ode(t+0.5*dt,   Y+0.5*dt*k2,    sysP, Om, Fw);
            k4 = ext_ode(t+dt,       Y+dt*k3,        sysP, Om, Fw);
            Y = Y + (dt/6)*(k1+2*k2+2*k3+k4);
            t = t + dt;
        end

        PhiT = reshape(Y(7:end),6,6);
        ev = eig(PhiT);
        max_mu(k) = max(abs(ev));
        is_stable(k) = max_mu(k) < tol_stable;
    end
end

function dY = ext_ode(t, Y, sysP, Omega, Fw)
% 状态+变分方程（严格对齐你之前的 Get_Floquet_Stability ext_ode）
    y = Y(1:6);
    Phi = reshape(Y(7:end),6,6);

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

    df12_ddx = (be1+al1) + 3*ga1*dx^2;
    df2g_dx2 = be2 + 3*ga2*x2^2;
    df2g_dv2 = 2*mu*ze2;

    A = zeros(6);
    A(1,2)=1;
    A(2,1)=-df12_ddx;  A(2,3)=+df12_ddx;                 A(2,6)=theta;
    A(3,4)=1;
    A(4,1)= df12_ddx/mu;
    A(4,3)= (-df12_ddx - df2g_dx2)/mu;
    A(4,4)= -df2g_dv2/mu;
    A(4,6)= -theta/mu;
    A(5,6)=1;

    if kap_e~=0
        A(6,2)=-theta/kap_e; A(6,4)=+theta/kap_e;
        A(6,5)=-kap_c/kap_e; A(6,6)=-sigma/kap_e;
    else
        s=max(abs(sigma),1e-12);
        A(6,2)=-50*theta/s;  A(6,4)=+50*theta/s;
        A(6,5)=-50*kap_c/s;  A(6,6)=-50;
    end

    dPhi = A*Phi;
    dY = [dy; dPhi(:)];
end

function [Fout,Aout,muout,stout] = bin_unique_by_F(F,A,mu,st,dF)
% F按 dF 分箱；每箱取“最保守点”：mu最大的那个
    F = F(:); A=A(:); mu=mu(:); st=st(:);
    b = round(F/dF);
    ub = unique(b,'stable');

    Fout = zeros(numel(ub),1);
    Aout = zeros(numel(ub),1);
    muout = zeros(numel(ub),1);
    stout = false(numel(ub),1);

    for i=1:numel(ub)
        idx = find(b==ub(i));
        [~,kmax] = max(mu(idx));
        k = idx(kmax);
        Fout(i)=F(k);
        Aout(i)=A(k);
        muout(i)=mu(k);
        stout(i)=st(k);
    end
end

function intervals = get_intervals_from_stability(F, is_stable, wantStable)
% 按F排序后，提取稳定/不稳定区间
% wantStable=false -> 不稳定区间
    [Fs, idx] = sort(F);
    st = is_stable(idx);
    target = (st == wantStable);

    intervals = [];
    if ~any(target), return; end

    k = 1;
    while k <= numel(Fs)
        if target(k)
            k0 = k;
            while k <= numel(Fs) && target(k)
                k = k + 1;
            end
            k1 = k - 1;
            intervals = [intervals; Fs(k0), Fs(k1)]; %#ok<AGROW>
        else
            k = k + 1;
        end
    end
end

function shade_intervals(ax, intervals)
% 用淡灰阴影标注区间
    if isempty(intervals), return; end
    yl = ylim(ax);
    for i=1:size(intervals,1)
        x1 = intervals(i,1); x2 = intervals(i,2);
        patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.85 0.85 0.85], 'EdgeColor','none', 'FaceAlpha',0.25);
    end
    ylim(ax, yl);
end

function [Fup, Fdown] = estimate_jump_from_hysteresis(Ff,Af,Fb,Ab)
% 用幅值最大突变估计 jump（够用、鲁棒）
    [Fup,~] = pick_jump(Ff,Af);
    [Fdown,~] = pick_jump(Fb,Ab);

    function [Fj, idx] = pick_jump(F,A)
        dA = abs(diff(A));
        [~,idx] = max(dA);
        Fj = F(idx+1);
    end
end
