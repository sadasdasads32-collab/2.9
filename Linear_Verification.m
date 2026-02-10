%% Linear Verification Script for Nondim_temp2 (Symmetric Coupling Version)
% 目的：验证 nondim_temp2 在物理一致性修正（sqrt(lam)）后是否正确
% 修正点：理论矩阵 C_lin 中的耦合项现在使用 theta = sqrt(lam) 对称分布

clear; clc; close all;

% --- 1. 参数设置 (System Parameters) ---
% 任意选取一组参数，但确保稳定
be1 = 1.0; 
be2 = 1.0; 
mu  = 0.5;
al1 = 0.01;   % Linear QZS stiffness part
ga1 = 0;      % 【关键】关闭非线性
ze1 = 0.05;
lam = 0.2;    % 物理机电耦合强度
kap_e = 0.05;
kap_c = 1.0;
sigma = 0.5;
ga2 = 0;      % 【关键】关闭非线性

% Pack parameters
sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2];

global Fw FixedOmega
Fw = 0.01;       % Excitation amplitude
FixedOmega = []; % Frequency sweep mode

% Frequency Range
Omega_List = 0.1:0.05:3.0;

% --- 2. 理论解计算 (Theoretical Linear FRF - Corrected) ---
% 计算对称耦合系数
theta = sqrt(lam); 

% System Matrices: M*q'' + C*q' + K*q = F
% Variables q = [x1; x2; Q]

M_lin = [1,    0,     0;
         0,    mu,    0;
         0,    0,     kap_e];

% 【核心修正】：阻尼/耦合矩阵 C_lin
% 机械方程1：... - theta*Q' ...  => C(1,3) = -theta
% 机械方程2：... + theta*Q' ...  => C(2,3) = +theta
% 电路方程3：... - theta*x1' + theta*x2' ... => C(3,1)=-theta, C(3,2)=+theta
C_lin = [0,          0,           -theta;
         0,          2*mu*ze1,    theta;
         -theta,     theta,       sigma]; 

K_lin = [ (be1+al1),    -(be1+al1),   0;
         -(be1+al1),    (be1+al1)+be2,0;
          0,            0,            kap_c];

F_vec = [Fw; 0; 0];

Amp_Theory = zeros(length(Omega_List), 3); % Store |x1|, |x2|, |Q|

for i = 1:length(Omega_List)
    w = Omega_List(i);
    % Dynamic Stiffness Matrix
    D = K_lin - w^2 * M_lin + 1i * w * C_lin;
    Response = D \ F_vec;
    Amp_Theory(i, :) = abs(Response)';
end

% --- 3. 数值解验证 (Numerical Validation via fsolve) ---
Amp_Num = zeros(length(Omega_List), 3);

% 启用解析雅可比加速
opts = optimoptions('fsolve','Display','off', ...
    'SpecifyObjectiveGradient',true, ...
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10);

fprintf('Starting Numerical Validation (Symmetric Check)...\n');
for i = 1:length(Omega_List)
    w = Omega_List(i);
    y0 = zeros(16,1);
    y0(16) = w; 
    
    % Solve using nondim_temp2
    solve_func = @(y) nondim_temp2([y(1:15); w], sysP);
    [y_sol, ~, exitflag] = fsolve(solve_func, zeros(15,1), opts);
    
    if exitflag > 0
        % Extract 1st harmonic Amplitude
        amp_x1 = norm(y_sol(2:3)); 
        amp_x2 = norm(y_sol(7:8));
        amp_q  = norm(y_sol(12:13));
        Amp_Num(i, :) = [amp_x1, amp_x2, amp_q];
    else
        Amp_Num(i, :) = [NaN, NaN, NaN];
    end
end

% --- 4. 绘图对比 ---
figure('Position',[100,100,1000,400]);
subplot(1,3,1); 
plot(Omega_List, Amp_Theory(:,1), 'k-', 'LineWidth', 2); hold on;
plot(Omega_List, Amp_Num(:,1), 'r--', 'LineWidth', 2);
title('X1 Amplitude'); xlabel('\Omega'); legend('Theory (Sym)','HBM-Num'); grid on;

subplot(1,3,2); 
plot(Omega_List, Amp_Theory(:,2), 'k-', 'LineWidth', 2); hold on;
plot(Omega_List, Amp_Num(:,2), 'r--', 'LineWidth', 2);
title('X2 Amplitude'); xlabel('\Omega'); grid on;

subplot(1,3,3); 
plot(Omega_List, Amp_Theory(:,3), 'k-', 'LineWidth', 2); hold on;
plot(Omega_List, Amp_Num(:,3), 'r--', 'LineWidth', 2);
title('Q Amplitude'); xlabel('\Omega'); grid on;

% Error Check
Max_Error = max(max(abs(Amp_Theory - Amp_Num)));
fprintf('Maximum Discrepancy: %e\n', Max_Error);

if Max_Error < 1e-6
    disp('✅ SUCCESS: Code matches SYMMETRIC linear theory perfectly.');
else
    disp('❌ FAILURE: Mismatch detected.');
end