%% Verify_Cubic_Truncation_Error_ExactVs3rd.m
% 三次截断误差验证（精确恢复力 vs 三阶近似）
% - 精确: F_h(z)=2k_h*(1 - L0/sqrt(L^2+z^2))*z
% - 三阶: F_h3(z)=2k_h*(1 - L0/L)*z + (k_h*L0/L^3)*z^3
% - 工作区间: 由 FRF 分支峰值位移确定（max|x12(t)|, max|x2(t)|）
% - 结果: 误差曲线 + 工作区间对齐的误差上界 e_max
%
% 依赖: FRF.m (你已上传), 以及你 FRF 所需的其他函数在路径下
% 注意: 本脚本只用于“恢复力截断误差”，不要求 nondim_temp2 内部改动

clear; clc; close all;

%% ---------------- 0) 基本设置 ----------------
Omega_min = 0.2;
Omega_max = 2.0;          % 先验证 0.2~2.0；若要全频段可改大
Nx = 3001;                % 误差曲线采样密度
Nt = 2048;                % 从 FRF 重构时域用于峰值估计

global Fw FixedOmega ParamMin ParamMax
FixedOmega = [];          % FRF: y(16)=Omega
Fw = 0.005;

% ====== 你的 QZS 参数（与你前面脚本一致） ======
mu   = 0.2;
beta = 2.0;
K1   = 1.0;
K2   = 0;
U    = 2.0;
L    = 4/9;

v = 2.5;
alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

% ====== nondim_temp2/FRF 所用参数包（与你保持一致） ======
P.be1 = 1.0;
P.al1 = alpha1 - P.be1;
P.be2 = alpha2;
P.ga1 = gamma1;
P.ga2 = gamma2;
P.mu  = mu;
P.ze1 = 0.05;

% 电路参数（不影响“单独的恢复力截断误差”，但影响 FRF 得到的峰值位移范围）
P.lam   = 0.18;
P.kap_e = 0.01;
P.kap_c = 0.01;
P.sigma = 0.43;

sysP = [P.be1; P.be2; P.mu; P.al1; P.ga1; P.ze1; P.lam; ...
        P.kap_e; P.kap_c; P.sigma; P.ga2];

fprintf('=== Cubic truncation error verification (Exact vs 3rd) ===\n');
fprintf('Omega range = [%.2f, %.2f], Fw=%.4g\n', Omega_min, Omega_max, Fw);

%% ---------------- 1) 跑 FRF 得到工作区间（峰值位移） ----------------
ParamMin = Omega_min - 0.05;
ParamMax = Omega_max;

x_frf = FRF(sysP);          % 16xN (含 Omega 行)
Om_all = x_frf(16,:).';

mask = (Om_all >= Omega_min) & (Om_all <= Omega_max);
x_frf = x_frf(:,mask);
Om = x_frf(16,:).';
N = numel(Om);

fprintf('FRF points used: N=%d\n', N);

x12_peak = zeros(N,1);
x2_peak  = zeros(N,1);

for i = 1:N
    W = Om(i);
    T = 2*pi/W;
    t = linspace(0,T,Nt).';

    x1c = x_frf(1:5,i);
    x2c = x_frf(6:10,i);

    x1 = recon_u_013(x1c,t,W);
    x2 = recon_u_013(x2c,t,W);
    x12 = x1 - x2;

    x12_peak(i) = max(abs(x12));
    x2_peak(i)  = max(abs(x2));
end

x12_peak_max = max(x12_peak);
x2_peak_max  = max(x2_peak);

fprintf('Operating range from FRF:\n');
fprintf('  max|x12(t)| = %.6g\n', x12_peak_max);
fprintf('  max|x2(t)|  = %.6g\n', x2_peak_max);

%% ---------------- 2) 构造精确恢复力与三阶近似 ----------------
% 你文中式(3)(4)：F_h(z;k_h,L,L0)
Fh_exact = @(z,kh,Ls,L0) 2*kh.*(1 - L0./sqrt(Ls.^2 + z.^2)).*z;
Fh_3rd   = @(z,kh,Ls,L0) 2*kh.*(1 - L0./Ls).*z + (kh.*L0./Ls.^3).*z.^3;

% 你文中总恢复力定义（含竖向线性弹簧并联）：
% 上层：Fr1(x)=k1*x + Fh(x;k3,L1,L01)
% 下层：Fr2(z2)=k2*z2 + Fh(z2;k4,L2,L02)
%
% 但你无量纲化后，用的是 (be1+al1)*x + ga1*x^3 这种形式
% 这里我们做“几何非线性截断误差”时，重点验证 Fh 的截断误差即可。
% 若你希望验证“总恢复力 Fr1/Fr2”的误差，也可以把 k1/k2 加进去（下面已给开关）。

include_linear_spring_in_total = true;   % =true: 比较 Fr1_exact vs Fr1_3rd（含线性弹簧）
                                         % =false: 只比较 Fh_exact vs Fh_3rd（纯几何项）

% ===== 关键：这里需要你把“几何参数映射到 (k_h,L,L0)” =====
% 你的无量纲写法里用了 K1,K2,L,U,v,beta 等。
% 若你想用“无量纲形式”直接比较，也行：令 k_h=K1, Ls=L, L0=(1-L) 或其它你论文定义。
%
% 你在 alpha_map 里出现 (1-L)/L，说明常见设定是 L0 = 1-L （无量纲长度）
% 这里按与你 alpha_map 一致的典型几何：Ls = L, L0 = 1-L 来做误差验证。
%

% 你论文的无量纲几何：L = (水平投影长度)/(原长)，因此 L0 = 1
L1  = L;     L01 = 1;     k3 = K1;   % 上层侧向弹簧
L2  = L;     L02 = 1;     k4 = K2;   % 下层侧向弹簧

% 竖向线性弹簧（无量纲）
k1v = v;        % 对应你 v = kv1/kv
k2v = beta;     % 对应 beta = kv2/kv

% 上层
Fr1_exact = @(x) (include_linear_spring_in_total)*k1v*x + Fh_exact(x,k3,L1,L01);
Fr1_3rd   = @(x) (include_linear_spring_in_total)*k1v*x + Fh_3rd  (x,k3,L1,L01);

% 下层
Fr2_exact = @(z) (include_linear_spring_in_total)*k2v*z + Fh_exact(z,k4,L2,L02);
Fr2_3rd   = @(z) (include_linear_spring_in_total)*k2v*z + Fh_3rd  (z,k4,L2,L02);

%% ---------------- 3) 在工作区间内计算误差曲线 ----------------
eps_den = 1e-14;

x12_grid = linspace(-x12_peak_max, x12_peak_max, Nx).';
x2_grid  = linspace(-x2_peak_max , x2_peak_max , Nx).';

Fex12 = Fr1_exact(x12_grid);
Fap12 = Fr1_3rd  (x12_grid);
err12 = abs(Fex12 - Fap12) ./ max(eps_den, abs(Fex12));  % relative

Fex2  = Fr2_exact(x2_grid);
Fap2  = Fr2_3rd  (x2_grid);
err2  = abs(Fex2 - Fap2) ./ max(eps_den, abs(Fex2));

err12_max = max(err12);
err2_max  = max(err2);

fprintf('\n=== Error bound within operating range ===\n');
fprintf('Upper (x12): e_max = %.6g (%.4f%%)\n', err12_max, err12_max*100);
fprintf('Lower (x2 ): e_max = %.6g (%.4f%%)\n', err2_max , err2_max*100);

%% ---------------- 4) 画图：恢复力对比 + 误差曲线（线性+log） ----------------
figure('Color','w','Position',[120,60,1100,820]);

% (a) 上层恢复力对比
subplot(2,2,1);
plot(x12_grid, Fex12, 'LineWidth',1.3); hold on;
plot(x12_grid, Fap12, '--', 'LineWidth',1.3);
grid on; box on;
xlabel('$x_{12}$','Interpreter','latex');
ylabel('$F(x)$','Interpreter','latex');
title('Upper restoring force: exact vs 3rd','Interpreter','latex');
legend({'$F_{\rm exact}$','$F_{3}$'},'Interpreter','latex','Location','best');

% (b) 上层误差（线性）
subplot(2,2,2);
plot(x12_grid, err12*100, 'LineWidth',1.3); hold on;
yline(err12_max*100,'--','LineWidth',1.2);
grid on; box on;
xlabel('$x_{12}$','Interpreter','latex');
ylabel('Relative error (\%)','Interpreter','latex');
title(sprintf('Upper error, $|x_{12}|\\le %.4g$, $e_{\\max}=%.3g\\%%$', ...
      x12_peak_max, err12_max*100),'Interpreter','latex');
% 自适应 y 轴（避免看成一条线）
ymax = max(err12)*100;
if ymax < 1e-3, ylim([0 1e-3]); else, ylim([0 1.1*ymax]); end

% (c) 下层恢复力对比
subplot(2,2,3);
plot(x2_grid, Fex2, 'LineWidth',1.3); hold on;
plot(x2_grid, Fap2, '--', 'LineWidth',1.3);
grid on; box on;
xlabel('$x_{2}$','Interpreter','latex');
ylabel('$F(x)$','Interpreter','latex');
title('Lower restoring force: exact vs 3rd','Interpreter','latex');
legend({'$F_{\rm exact}$','$F_{3}$'},'Interpreter','latex','Location','best');

% (d) 下层误差（半对数更清楚）
subplot(2,2,4);
semilogy(x2_grid, max(err2*100,1e-12), 'LineWidth',1.3); hold on;
yline(max(err2_max*100,1e-12),'--','LineWidth',1.2);
grid on; box on;
xlabel('$x_{2}$','Interpreter','latex');
ylabel('Relative error (\%), log scale','Interpreter','latex');
title(sprintf('Lower error, $|x_2|\\le %.4g$, $e_{\\max}=%.3g\\%%$', ...
      x2_peak_max, err2_max*100),'Interpreter','latex');
legend({'$e(x)$','$e_{\max}$'},'Interpreter','latex','Location','best');

%% ---------------- 5) 误差上界与工作区间对齐展示（Ω 维度） ----------------
figure('Color','w','Position',[1250,140,900,520]);

yyaxis left
plot(Om, x12_peak, 'LineWidth',1.3); hold on;
plot(Om, x2_peak , 'LineWidth',1.3);
ylabel('Peak displacement','Interpreter','latex');

yyaxis right
plot(Om, err12_max*100*ones(size(Om)), '--', 'LineWidth',1.3);
ylabel('$e_{\max}$ in operating range (\%)','Interpreter','latex');

grid on; box on;
xlabel('$\Omega$','Interpreter','latex');
title('Operating range aligned with truncation-error upper bound','Interpreter','latex');
legend({'$\max|x_{12}(t)|$','$\max|x_2(t)|$','$e_{\max}$ (upper)'}, ...
       'Interpreter','latex','Location','best');

%% =================== helper: HB 0/1/3 reconstruction ===================
function u = recon_u_013(cfs, t, W)
    u = cfs(1) ...
      + cfs(2)*cos(W*t) + cfs(3)*sin(W*t) ...
      + cfs(4)*cos(3*W*t) + cfs(5)*sin(3*W*t);
end