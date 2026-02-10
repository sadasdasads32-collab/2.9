%% Linear Verification Script for Nondim_temp2
% 目的：验证 nondim_temp2 在 ga1=0, ga2=0 时是否与线性理论解一致
% 拓扑：Standard Two-Stage Raft (T1)

clear; clc; close all;

% --- 1. 参数设置 (System Parameters) ---
% 任意选取一组参数，但确保稳定
be1 = 1.0; 
be2 = 1.0; 
mu  = 0.5;
al1 = 0.01;   % Linear QZS stiffness part
ga1 = 0;      % 【关键】关闭非线性
ze1 = 0.05;
lam = 0.2;
kap_e = 0.05;
kap_c = 1.0;
sigma = 0.5;
ga2 = 0;      % 【关键】关闭非线性

% Pack parameters for nondim_temp2
% [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2];

global Fw FixedOmega
Fw = 0.01;       % Excitation amplitude
FixedOmega = []; % Frequency sweep mode

% Frequency Range
Omega_List = 0.1:0.05:3.0;

% --- 2. 理论解计算 (Theoretical Linear FRF) ---
% System Matrix form: M*q'' + C*q' + K*q = F
% Variables q = [x1; x2; Q]
% Equations:
% 1) x1'' + (be1+al1)(x1-x2) - lam*Q' = f
% 2) mu*x2'' + 2*mu*ze1*x2' + be2*x2 + (be1+al1)(x2-x1) + lam*Q' = 0
% 3) kap_e*Q'' + sigma*Q' + kap_c*Q - (x1' - x2') = 0

M_lin = [1,    0,     0;
         0,    mu,    0;
         0,    0,     kap_e];

C_lin = [0,          0,           -lam;
         0,          2*mu*ze1,    lam;
         -1,         1,           sigma]; % Note: Eq3 term -(x1'-x2') -> -x1' + x2'

K_lin = [ (be1+al1),    -(be1+al1),   0;
         -(be1+al1),    (be1+al1)+be2,0;
          0,            0,            kap_c];

F_vec = [Fw; 0; 0];

Amp_Theory = zeros(length(Omega_List), 3); % Store |x1|, |x2|, |Q|

for i = 1:length(Omega_List)
    w = Omega_List(i);
    % Dynamic Stiffness Matrix: D = K - w^2*M + j*w*C
    D = K_lin - w^2 * M_lin + 1i * w * C_lin;
    
    % Linear Solve
    Response = D \ F_vec;
    Amp_Theory(i, :) = abs(Response)';
end

% --- 3. 数值解验证 (Numerical Validation via fsolve) ---
Amp_Num = zeros(length(Omega_List), 3);
y_guess = zeros(16,1); 

% Optimization options: Enable Jacobian
opts = optimoptions('fsolve','Display','off', ...
    'SpecifyObjectiveGradient',true, ... % <--- 【关键】启用解析雅可比
    'FunctionTolerance',1e-10, 'StepTolerance',1e-10);

fprintf('Starting Numerical Validation...\n');
for i = 1:length(Omega_List)
    w = Omega_List(i);
    
    % Update guess with linear theory (perfect warm start)
    % Just to check if residual is zero at exact solution
    % Or solve from scratch to test convergence
    
    % Construct guess from theory (to be robust)
    % Theory: X = A - jB.  HBM: cos coeff = A, sin coeff = B
    % Response is X * e^{jwt} = (A+jB)*(cos+jsin) = Acos - Bsin + j(...)
    % My code usually assumes x = u(2)*cos + u(3)*sin
    % Phasor X = u(2) - j*u(3).  => u(2)=real(X), u(3)=-imag(X)
    
    % Let's run fsolve from zero to test robustness
    y0 = zeros(16,1);
    y0(16) = w; % Current frequency
    
    % Solve
    % Pass w as parameter inside y (continuation param style)
    solve_func = @(y) nondim_temp2([y(1:15); w], sysP);
    
    [y_sol, ~, exitflag] = fsolve(solve_func, zeros(15,1), opts);
    
    if exitflag > 0
        % Extract 1st harmonic Amplitude: sqrt(u2^2 + u3^2)
        % x1 is indices 1:5. 1st harmonic is 2,3
        amp_x1 = norm(y_sol(2:3)) / 1; % HBM defines x = u2*c + u3*s. Amp = sqrt(u2^2+u3^2)
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
title('X1 Amplitude'); xlabel('\Omega'); legend('Theory','HBM-Num'); grid on;

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
fprintf('Maximum Discrepancy between Theory and Code: %e\n', Max_Error);

if Max_Error < 1e-6
    disp('✅ SUCCESS: Code matches linear theory perfectly.');
else
    disp('❌ FAILURE: Significant mismatch detected. Check equations again.');
end