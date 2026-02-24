%% Run_Step1c_Theory_Damped_vs_FRF_WindowPeaks_FIXED.m


clear; clc; close all;

%% -----------------------------
% 1) 你的参数（保持不动）
% -----------------------------
P.be1 = 1;
P.be2 = 1.0;
P.mu  = 0.2;
P.al1 = 0;

P.ga1 = 0.0;
P.ga2 = 0.0;

% 你的数值模型：下层对地阻尼=ze1；lam=0时 nondim_temp2 内部启用 zeta12=0.020412
P.ze1   = 0.11183;
P.lam   = 0.0;
P.kap_e = 0.0;
P.kap_c = 0.0;
P.sigma = 0.0;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.005;
FixedOmega = [];

%% -----------------------------
% 2) 理论：工程阻尼上下都 0.05 -> 无量纲 zeta2, zeta12（只用于理论）
% -----------------------------
zeta2_eng  = 0.05;
zeta12_eng = 0.05;

mu = P.mu;
b2 = P.be2;
a1 = P.al1;

zeta2  = zeta2_eng  * sqrt(b2/mu);
zeta12 = zeta12_eng * sqrt(mu*(1+a1)/(1+mu));

fprintf('================================================\n');
fprintf('THEORY (eng damp up/down = 0.05) vs NUMERICAL FRF(sysP)\n');
fprintf('Converted nondim theory: zeta2=%.6f, zeta12=%.6f\n', zeta2, zeta12);
fprintf('Your numerical inputs: ze1=%.6f, lam=%.6f\n', P.ze1, P.lam);
fprintf('================================================\n\n');

%% -----------------------------
% 3) 线性 M, K
% -----------------------------
k12 = (P.be1 + P.al1);
k2  = P.be2;

M = diag([1, P.mu]);
K = [ k12,    -k12;
     -k12, k12 + k2];

wn = sqrt(eig(M\K));
wn = sort(real(wn));

fprintf('--- Undamped wn (from M\\K) ---\n');
fprintf('  wn1 = %.6f\n', wn(1));
fprintf('  wn2 = %.6f\n', wn(2));
fprintf('--------------------------------\n\n');

%% -----------------------------
% 4) 理论阻尼矩阵 C（相对阻尼 + 对地阻尼）
% -----------------------------
C = [ 2*zeta12,              -2*zeta12;
     -2*zeta12,  2*zeta12 + 2*P.mu*zeta2 ];

%% -----------------------------
% 5) 含阻尼特征频率 wd
% -----------------------------
Z = zeros(2); I = eye(2);
A = [ Z,      I;
     -M\K,  -M\C ];

lamA = eig(A);
lamA = lamA(imag(lamA) > 1e-9);
wd = sort(imag(lamA));
wd = wd(1:min(2,end));

fprintf('--- Damped modal frequencies wd = Im(eig(A)) ---\n');
fprintf('  wd1 = %.6f\n', wd(1));
fprintf('  wd2 = %.6f\n', wd(2));
fprintf('-----------------------------------------------\n\n');

%% -----------------------------
% 6) 理论解析 FRF（动态刚度法）
% -----------------------------
Omega_th = linspace(0.01, max(5, 1.2*wn(2)), 12000);  % 再加密一点
Fvec = [Fw; 0];

X1_th = zeros(size(Omega_th));
X2_th = zeros(size(Omega_th));

for k = 1:numel(Omega_th)
    Om = Omega_th(k);
    D = -Om^2*M + 1i*Om*C + K;
    X = D \ Fvec;
    X1_th(k) = X(1);
    X2_th(k) = X(2);
end

eps_db = 1e-14;
Amp1_th_dB = 20*log10(abs(X1_th) + eps_db);
Amp2_th_dB = 20*log10(abs(X2_th) + eps_db);

%% -----------------------------
% 7) 数值 FRF（你的 HBM+延拓）
% -----------------------------
x_num = FRF(sysP);

Omega_num_raw = x_num(16,:);
Amp1_num_raw  = sqrt(x_num(2,:).^2 + x_num(3,:).^2);
Amp2_num_raw  = sqrt(x_num(7,:).^2 + x_num(8,:).^2);

Amp1_num_dB_raw = 20*log10(Amp1_num_raw + eps_db);
Amp2_num_dB_raw = 20*log10(Amp2_num_raw + eps_db);

% ===== 关键修正：按 Omega 排序（延拓点序可能乱）=====
[Omega_num, idxSort] = sort(Omega_num_raw, 'ascend');
Amp1_num_dB = Amp1_num_dB_raw(idxSort);
Amp2_num_dB = Amp2_num_dB_raw(idxSort);

%% -----------------------------
% 8) 绘图（只画曲线，不标记峰值点）
% -----------------------------
figure('Color','w','Position',[120 120 1200 450]);

subplot(1,2,1);
plot(Omega_th,  Amp1_th_dB,  'k-',  'LineWidth', 2); hold on;
plot(Omega_num, Amp1_num_dB, 'b--', 'LineWidth', 1.5);
grid on; axis tight; yline(0,'-');
xlabel('\Omega'); ylabel('|X_1| (dB)');
title('X1: Theory vs Numerical');
legend('Theory','Numerical','Location','best');

subplot(1,2,2);
plot(Omega_th,  Amp2_th_dB,  'k-',  'LineWidth', 2); hold on;
plot(Omega_num, Amp2_num_dB, 'r--', 'LineWidth', 1.5);
grid on; axis tight; yline(0,'-');
xlabel('\Omega'); ylabel('|X_2| (dB)');
title('X2: Theory vs Numerical');
legend('Theory','Numerical','Location','best');