%% ODE_MultiIC_Convergence.m
% 固定 sysP, Omega_test, Fw：多初值 ODE45 收敛性/多稳态检测
clc; close all;

global Fw
Fw = 0.005;

% ======== 填你的 sysP ========
P.be1=1.0; P.be2=0.1; P.mu=0.2; P.al1=-0.95;
k1=1; k2=0.8; L=4/9; U=2;
P.ga1=k1/(U^2*L^3);
P.ga2=k2/(U^2*L^3);
P.ze1=0.05;
P.lam=0.18; P.kap_e=0.2; P.kap_c=0.5; P.sigma=1.0;

sysP = [P.be1; P.be2; P.mu; P.al1; P.ga1; P.ze1; ...
        P.lam; P.kap_e; P.kap_c; P.sigma; P.ga2];

Omega_test = 0.80;
sgn_em = +1;

% ======== 多初值设置 ========
N_ic = 30;                 % 初值个数（建议 30~100）
scale_x = 5e-3;            % 位移初值尺度
scale_v = 5e-3;            % 速度初值尺度
scale_q = 5e-3;            % 电荷初值尺度
scale_qd= 5e-3;            % 电荷导数初值尺度

T_periods = 350;
N_last = 30;
T1 = 2*pi/Omega_test;
tspan = [0, T_periods*T1];

opts = odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',T1/50);

A1x_all = nan(N_ic,1);
ok_all  = false(N_ic,1);

fprintf('==== Multi-IC ODE45 @ Omega=%.4f, Fw=%.4f ====\n', Omega_test, Fw);

% 存几个示例时域用于画图
storeK = min(8, N_ic);
t_store = cell(storeK,1);
x_store = cell(storeK,1);

for i=1:N_ic
    % 随机初值
    y0 = zeros(6,1);
    y0(1) = scale_x * (2*rand-1);
    y0(2) = scale_v * (2*rand-1);
    y0(3) = scale_x * (2*rand-1);
    y0(4) = scale_v * (2*rand-1);
    y0(5) = scale_q * (2*rand-1);
    y0(6) = scale_qd* (2*rand-1);

    try
        [t,y] = ode45(@(t,y) odesys_nd(t,y,sysP,Omega_test,Fw,sgn_em), tspan, y0, opts);

        I = find(t > t(end) - N_last*T1);
        t_s = t(I); x1_s = y(I,1);
        x = x1_s - mean(x1_s);

        C = 2/numel(t_s) * sum(x .* cos(Omega_test*t_s));
        S = 2/numel(t_s) * sum(x .* sin(Omega_test*t_s));
        A1x = hypot(C,S);

        A1x_all(i) = A1x;
        ok_all(i) = isfinite(A1x);

        if i<=storeK
            t_store{i} = t_s;
            x_store{i} = x1_s;
        end

        fprintf('IC %2d/%2d OK: A1x=%.6g\n', i, N_ic, A1x);

    catch ME
        fprintf('IC %2d/%2d FAIL: %s\n', i, N_ic, ME.message);
        ok_all(i)=false;
    end
end

A = A1x_all(ok_all);
fprintf('\n==== Summary ====\n');
fprintf('Success: %d/%d\n', nnz(ok_all), N_ic);
fprintf('A1x mean=%.6g, std=%.3g\n', mean(A), std(A));

% 简单聚类：按相对阈值分组（审稿够用）
tol = max(1e-6, 1e-3*mean(A));  % 0.1% 的幅值尺度
A_sorted = sort(A);
clusters = {};
if ~isempty(A_sorted)
    cur = A_sorted(1);
    grp = cur;
    for k=2:numel(A_sorted)
        if abs(A_sorted(k)-cur) <= tol
            grp(end+1)=A_sorted(k);
        else
            clusters{end+1}=grp; %#ok<AGROW>
            grp = A_sorted(k);
        end
        cur = A_sorted(k);
    end
    clusters{end+1}=grp;
end

fprintf('Cluster count = %d (tol=%.2e)\n', numel(clusters), tol);
for c=1:numel(clusters)
    fprintf('  Cluster %d: count=%d, mean=%.6g\n', c, numel(clusters{c}), mean(clusters{c}));
end

% ======== 绘图 ========
figure('Color','w','Name','A1x scatter across ICs');
plot(find(ok_all), A1x_all(ok_all), 'o'); grid on;
xlabel('IC index'); ylabel('A_{1x}');
title(sprintf('ODE Multi-IC steady amplitude @ \\Omega=%.4f (tol=%.2e)', Omega_test, tol));

figure('Color','w','Name','A1x histogram');
histogram(A, 12); grid on;
xlabel('A_{1x}'); ylabel('count');
title('Histogram of steady-state A_{1x}');

figure('Color','w','Name','Some steady-state time traces'); hold on; grid on;
for i=1:storeK
    if ~isempty(t_store{i})
        plot(t_store{i}, x_store{i});
    end
end
xlabel('t'); ylabel('x_1');
title(sprintf('Last %d cycles time traces (first %d ICs stored)', N_last, storeK));

%% ==== ODE 方程（同你验证脚本）====
function dydt = odesys_nd(t, y, sysP, Omega, Fw, sgn_em)
    x1=y(1); v1=y(2); x2=y(3); v2=y(4); q=y(5); qd=y(6);
    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    dx  = x1-x2;
    f12 = (be1+al1)*dx + ga1*dx^3;
    f2g = be2*x2 + ga2*x2^3 + 2*mu*ze2*v2;

    th = sqrt(max(lam,0));

    x1dd = -f12 + sgn_em*(th*qd) + Fw*cos(Omega*t);
    x2dd = ( f12 - f2g - sgn_em*(th*qd) ) / mu;

    if kap_e == 0
        tiny = 1e-12;
        qd_new = (-kap_c*q - th*(v1 - v2)) / max(abs(sigma), tiny);
        qdd = (qd_new - qd)*50;
    else
        qdd = (-sigma*qd - kap_c*q - th*(v1 - v2)) / kap_e;
    end

    dydt = [v1; x1dd; v2; x2dd; qd; qdd];
end
