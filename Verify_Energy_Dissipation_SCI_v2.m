%% Verify_Energy_Dissipation_SCI_v2.m
% Strong rebuttal-grade verification:
% (1) Cycle-averaged power balance: <P_in> ≈ <P_mech> + <P_elec>
% (2) Circuit source consistency: e_LHS = kap_e*q'' + sigma*q' + kap_c*q  ≈  theta*(x1'-x2')
% (3) EM force consistency from mech eq: F_em_LHS ≈ theta*q'
%
% Dependencies: FRF.m, nondim_temp2.m (your final consistent version)

clear; clc; close all;

%% -------- Case setting --------
Current_Omega = 0.70;
Current_Lam   = 0.40;

% Parameters (match your sweep)
P.be1   = 1.0;
P.be2   = 1.0;
P.mu    = 0.2;
P.al1   = 0.0;
P.ga1   = 0.1;
P.ze1   = 0.05;
P.lam   = Current_Lam;
P.kap_e = 0.05;
P.kap_c = 1.0;
P.sigma = 2.0;
P.ga2   = 0.2;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
        P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.05;
FixedOmega = [];   % FRF sweep mode for initial guess

theta = sqrt(max(P.lam,0));

fprintf('================================================\n');
fprintf('ENERGY VERIFICATION (SCI v2): lambda=%.2f, Omega=%.3f\n', P.lam, Current_Omega);
fprintf('================================================\n');

%% -------- 1) initial guess from FRF --------
x_res = FRF(sysP);
if isempty(x_res) || size(x_res,1) < 16
    error('FRF output invalid. Check FRF.m and dependencies.');
end
Omega_grid = x_res(16,:);
[~, idx0] = min(abs(Omega_grid - Current_Omega));
y_init15 = x_res(1:15, idx0);

fprintf('Initial guess from FRF at Omega=%.4f\n', Omega_grid(idx0));

%% -------- 2) solve fixed Omega HBM --------
solve_func = @(y15) nondim_temp2([y15(:); Current_Omega], sysP);

opt = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxIterations',500, ...
    'MaxFunctionEvaluations',20000);

[y_sol15, fval, exitflag] = fsolve(solve_func, y_init15, opt);
fprintf('fsolve exitflag = %d, ||res||_inf = %.3e\n', exitflag, norm(fval,inf));

%% -------- 3) reconstruct time histories --------
Omega = Current_Omega;
T = 2*pi/Omega;
Nt = 800;
tau = linspace(0, T, Nt);

x1cfs = y_sol15(1:5);
x2cfs = y_sol15(6:10);
qcfs  = y_sol15(11:15);

% signals
x1  = recon_u (x1cfs, tau, Omega);
x2  = recon_u (x2cfs, tau, Omega);
q   = recon_u (qcfs,  tau, Omega);

x1p = recon_up(x1cfs, tau, Omega);
x2p = recon_up(x2cfs, tau, Omega);
qp  = recon_up(qcfs,  tau, Omega);          % qp = Q' = current (nondim)

x1pp = recon_upp(x1cfs, tau, Omega);
x2pp = recon_upp(x2cfs, tau, Omega);
qpp  = recon_upp(qcfs,  tau, Omega);

Fexc = Fw*cos(Omega*tau);

% relative
x12   = x1 - x2;
x12p  = x1p - x2p;

%% -------- 4) power terms --------
Pin  = Fexc .* x1p;

% mechanical damping on raft: 2*mu*zeta2 * x2'^2
c2 = 2*P.mu*P.ze1;
Pmech = c2*(x2p.^2);

% resistor dissipation
Pelec = P.sigma*(qp.^2);

% cycle averages
Mean_Pin   = trapz(tau, Pin)/T;
Mean_Pmech = trapz(tau, Pmech)/T;
Mean_Pelec = trapz(tau, Pelec)/T;

Bal_err = abs(Mean_Pin - (Mean_Pmech+Mean_Pelec)) / max(1e-12,abs(Mean_Pin)) * 100;

%% -------- 5) HARD checks (non-circular) --------
% (A) Circuit source voltage consistency:
% From circuit eq: kap_e*q'' + sigma*q' + kap_c*q = e
e_LHS = P.kap_e*qpp + P.sigma*qp + P.kap_c*q;
e_RHS = theta * x12p;
e_err_rms = rms(e_LHS - e_RHS) / max(1e-12, rms(e_RHS)) * 100;  % %

% (B) EM force consistency from mechanical eq (use R1 structure):
% Your nondim R1 is: x1'' + (be1+al1)*x12 + ga1*x12^3 + force_em = Fexc
% => force_em_LHS = Fexc - [x1'' + (be1+al1)*x12 + ga1*x12^3]
% => should match force_em_RHS = theta*qp
force_em_LHS = Fexc - (x1pp + (P.be1+P.al1)*x12 + P.ga1*(x12.^3));
force_em_RHS = theta * qp;
fem_err_rms = rms(force_em_LHS - force_em_RHS) / max(1e-12, rms(force_em_RHS)) * 100;

%% -------- 6) Report --------
fprintf('\n--- Cycle-averaged Power Balance ---\n');
fprintf('<P_in>           = %.8e\n', Mean_Pin);
fprintf('<P_mech_damp>    = %.8e\n', Mean_Pmech);
fprintf('<P_elec_res>     = %.8e\n', Mean_Pelec);
fprintf('<P_mech+P_elec>  = %.8e\n', Mean_Pmech+Mean_Pelec);
fprintf('Balance error    = %.6f %%\n', Bal_err);

fprintf('\n--- Non-circular consistency checks ---\n');
fprintf('Circuit voltage check:  e_LHS vs e_RHS, RMS error = %.6f %%\n', e_err_rms);
fprintf('EM force check:         fem_LHS vs fem_RHS, RMS error = %.6f %%\n', fem_err_rms);

pass = (Bal_err < 0.5) && (e_err_rms < 1.0) && (fem_err_rms < 1.0) && (Mean_Pelec > 0) && (P.sigma > 0);
if pass
    fprintf('\n[PASS] Strong evidence: model is reciprocal + resistor dissipates energy + power balance closes.\n');
else
    fprintf('\n[WARN] One or more checks failed. Re-check branch / parameters / sign conventions.\n');
end

%% -------- 7) Plots --------
figure('Color','w','Position',[160,80,980,860]);

subplot(3,1,1);
plot(tau, x1, 'b', 'LineWidth',1.3); hold on;
plot(tau, x2, 'r--', 'LineWidth',1.3);
grid on; xlim([0,T]);
ylabel('x');
legend('x_1','x_2');
title(sprintf('Time histories (\\Omega=%.3f, \\lambda=%.2f)', Omega, P.lam));

subplot(3,1,2);
plot(tau, Pin, 'k', 'LineWidth',1.1); hold on;
plot(tau, Pmech, 'r-.', 'LineWidth',1.1);
plot(tau, Pelec, 'g', 'LineWidth',1.4);
grid on; xlim([0,T]);
ylabel('Power');
legend('P_{in}','P_{mech}','P_{elec}');
title(sprintf('Power terms (balance err = %.4g%%)', Bal_err));

subplot(3,1,3);
plot(tau, e_LHS, 'm', 'LineWidth',1.1); hold on;
plot(tau, e_RHS, 'c--', 'LineWidth',1.1);
grid on; xlim([0,T]);
xlabel('\tau'); ylabel('e');
legend('e_{LHS}=\\kappa_eQ''''+\\sigmaQ''+\\kappa_cQ','e_{RHS}=\\theta(\\xi_1''-\\xi_2'')','Location','best');
title(sprintf('Circuit source consistency (RMS err = %.4g%%)', e_err_rms));

figure('Color','w','Position',[1180,140,880,420]);
plot(tau, force_em_LHS, 'k', 'LineWidth',1.1); hold on;
plot(tau, force_em_RHS, 'g--', 'LineWidth',1.1);
grid on; xlim([0,T]);
xlabel('\tau'); ylabel('F_{em}');
legend('F_{em,LHS} from mech eq','F_{em,RHS}=\\thetaQ''','Location','best');
title(sprintf('EM force consistency (RMS err = %.4g%%)', fem_err_rms));

%% -------- helper funcs --------
function u = recon_u(cfs, t, W)
u = cfs(1) + cfs(2)*cos(W*t) + cfs(3)*sin(W*t) + cfs(4)*cos(3*W*t) + cfs(5)*sin(3*W*t);
end

function up = recon_up(cfs, t, W)
up = cfs(2)*(-W*sin(W*t)) + cfs(3)*(W*cos(W*t)) + cfs(4)*(-3*W*sin(3*W*t)) + cfs(5)*(3*W*cos(3*W*t));
end

function upp = recon_upp(cfs, t, W)
% second derivative wrt nondim time
upp = cfs(2)*(-W^2*cos(W*t)) + cfs(3)*(-W^2*sin(W*t)) + cfs(4)*(-(3*W)^2*cos(3*W*t)) + cfs(5)*(-(3*W)^2*sin(3*W*t));
end
