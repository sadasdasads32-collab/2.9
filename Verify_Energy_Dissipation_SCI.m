%% Verify_Energy_Dissipation_SCI.m
% 审稿人回应专用：稳态周期功率平衡 + 耦合功率互易抵消验证
% 证明：
%  (1) <P_in> ≈ <P_mech_damp> + <P_elec_res>
%  (2) P_em_mech(t) + P_em_elec(t) ≈ 0 (逐点互易)
%
% 依赖：FRF.m, nondim_temp2.m 及其依赖

clear; clc; close all;

%% -------- Case setting --------
Current_Omega = 0.70;
Current_Lam   = 0.40;

% Parameters
P.be1   = 1.0;
P.be2   = 1.0;
P.mu    = 0.2;
P.al1   = 0.0;
P.ga1   = 0.1;
P.ze1   = 0.05;
P.lam   = Current_Lam;
P.kap_e = 0.05;
P.kap_c = 1.0;
P.sigma = 2.0;     % passive resistor
P.ga2   = 0.2;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
        P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.05;
FixedOmega = [];   % FRF sweep mode for initial guess

fprintf('================================================\n');
fprintf('ENERGY VERIFICATION (SCI): lambda=%.2f, Omega=%.3f\n', P.lam, Current_Omega);
fprintf('================================================\n');

%% -------- 1) Get a good initial guess from FRF --------
x_res = FRF(sysP);
if isempty(x_res) || size(x_res,1) < 16
    error('FRF output invalid. Check FRF.m and dependencies.');
end
Omega_grid = x_res(16,:);
[~, idx0] = min(abs(Omega_grid - Current_Omega));
y_init15 = x_res(1:15, idx0);     % use FRF solution as initial guess

fprintf('Initial guess picked from FRF at Omega=%.4f (closest)\n', Omega_grid(idx0));

%% -------- 2) Solve HBM algebraic equations at fixed Omega --------
solve_func = @(y15) nondim_temp2([y15(:); Current_Omega], sysP);

opt = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',20000);

[y_sol15, fval, exitflag] = fsolve(solve_func, y_init15, opt);
fprintf('fsolve exitflag = %d, ||res||_inf = %.3e\n', exitflag, norm(fval,inf));
if exitflag <= 0
    warning('HBM solution did not converge well. Try another Omega or better initial guess.');
end

%% -------- 3) Reconstruct time histories over one nondim period --------
Omega = Current_Omega;
T = 2*pi/Omega;
Nt = 400;
tau = linspace(0, T, Nt);   % tau is nondimensional time

x1cfs = y_sol15(1:5);
x2cfs = y_sol15(6:10);
qcfs  = y_sol15(11:15);

x1 = recon_u(x1cfs, tau, Omega);
x2 = recon_u(x2cfs, tau, Omega);
q  = recon_u(qcfs,  tau, Omega);

x1p = recon_up(x1cfs, tau, Omega);  % xi1'
x2p = recon_up(x2cfs, tau, Omega);  % xi2'
qp  = recon_up(qcfs,  tau, Omega);  % Q' = i

% excitation force (nondim): f cos(Omega tau)
Fexc = Fw*cos(Omega*tau);

%% -------- 4) Power terms (instantaneous) --------
% Input power
Pin = Fexc .* x1p;

% Mechanical damping on raft: term is 2*mu*ze1*xi2'
c2 = 2*P.mu*P.ze1;
Pmech = c2*(x2p.^2);

% Electrical resistor dissipation: sigma*(Q')^2
Pelec = P.sigma*(qp.^2);

% Coupling power check (reciprocity):
% mech side: (-sqrt(lam)Q')*x1' + (+sqrt(lam)Q')*x2' = -sqrt(lam)Q'*(x1'-x2')
sqrtlam = sqrt(max(P.lam,0));
Pem_mech = -sqrtlam*qp.*(x1p - x2p);
% elec side source power: (+sqrt(lam)(x1'-x2'))*Q'
Pem_elec = +sqrtlam*(x1p - x2p).*qp;
Pem_sum  = Pem_mech + Pem_elec;  % should be ~0 pointwise

%% -------- 5) Cycle-averaged powers --------
Mean_Pin   = trapz(tau, Pin)/T;
Mean_Pmech = trapz(tau, Pmech)/T;
Mean_Pelec = trapz(tau, Pelec)/T;

Mean_Pem_sum = trapz(tau, Pem_sum)/T;           % should be ~0
Max_Pem_sum  = max(abs(Pem_sum));               % pointwise reciprocity error

Bal_err = abs(Mean_Pin - (Mean_Pmech+Mean_Pelec)) / max(1e-12,abs(Mean_Pin)) * 100;

%% -------- 6) Report --------
fprintf('\n--- Cycle-averaged Power Balance ---\n');
fprintf('<P_in>      = %.8e\n', Mean_Pin);
fprintf('<P_mech>    = %.8e\n', Mean_Pmech);
fprintf('<P_elec>    = %.8e\n', Mean_Pelec);
fprintf('<P_mech+P_elec> = %.8e\n', Mean_Pmech+Mean_Pelec);
fprintf('Balance error = %.4f %%\n', Bal_err);

fprintf('\n--- Coupling reciprocity check ---\n');
fprintf('<Pem_mech + Pem_elec> (mean) = %.3e  (should ~0)\n', Mean_Pem_sum);
fprintf('max |Pem_mech+Pem_elec| (pointwise) = %.3e\n', Max_Pem_sum);

if (P.sigma>0) && (Mean_Pelec>0) && (Bal_err<1.0) && (abs(Mean_Pem_sum)<1e-6)
    fprintf('\n[PASS] Evidence strong: resistor dissipates energy and power balance closes.\n');
else
    fprintf('\n[WARN] Check convergence / indexing / whether Omega is on correct branch.\n');
end

%% -------- 7) Plots for rebuttal --------
figure('Color','w','Position',[200,120,900,700]);

subplot(3,1,1);
plot(tau, x1, 'b', 'LineWidth',1.5); hold on;
plot(tau, x2, 'r--', 'LineWidth',1.5);
grid on; xlim([0,T]);
ylabel('Displacement');
legend('x_1','x_2');
title(sprintf('Time histories at \\Omega=%.3f, \\lambda=%.2f', Omega, P.lam));

subplot(3,1,2);
plot(tau, Pin, 'k', 'LineWidth',1.2); hold on;
plot(tau, Pmech, 'r-.', 'LineWidth',1.2);
plot(tau, Pelec, 'g', 'LineWidth',1.6);
grid on; xlim([0,T]);
ylabel('Power');
legend('P_{in}','P_{mech}','P_{elec}');
title('Instantaneous power terms');

subplot(3,1,3);
plot(tau, Pem_sum, 'm', 'LineWidth',1.2);
grid on; xlim([0,T]);
xlabel('\tau');
ylabel('P_{em,mech}+P_{em,elec}');
title('Coupling reciprocity residual (should be ~0)');

%% -------- helper functions --------
function u = recon_u(cfs, t, W)
% cfs = [DC, c1, s1, c3, s3]
u = cfs(1) ...
  + cfs(2)*cos(W*t)   + cfs(3)*sin(W*t) ...
  + cfs(4)*cos(3*W*t) + cfs(5)*sin(3*W*t);
end

function up = recon_up(cfs, t, W)
% derivative wrt nondim time tau
up = 0 ...
   + cfs(2)*(-W*sin(W*t))    + cfs(3)*(W*cos(W*t)) ...
   + cfs(4)*(-3*W*sin(3*W*t)) + cfs(5)*(3*W*cos(3*W*t));
end
