%% Verify_Energy_Consistency_Current.m
% Strong verification (HBM-consistent, rebuttal-grade)
%  (1) Circuit equation residual check in time domain (R3(t) RMS)
%  (2) EM force consistency: F_em from mech eq vs theta*q'
%  (3) Cycle-averaged power balance: <Pin> ≈ <Pmech_damp> + <Pres>
%
% IMPORTANT SIGN CONVENTION (CORRECTED):
%   Your updated nondim_temp2 uses the physically correct Lenz's law:
%     R3 = (kap_e*q'' + sigma*q' + kap_c*q) - theta*(x1' - x2') = 0
%   => (kap_e*q'' + sigma*q' + kap_c*q) = + theta*(x1' - x2')
%   The checks below now follow THIS physically correct convention.

clear; clc; close all;

%% ---------- User case ----------
Omega_test = 0.70;      % choose the point you want to verify

global Fw FixedOmega
Fw = 0.005;
FixedOmega = [];        % keep FRF in sweep mode for initial guess

% ---- Your CURRENT parameter set example (replace as needed) ----
k1=1; k2=0.8;
L=4/9; U=2;

P.be1 = 1.0;
P.mu  = 0.2;
P.be2 = 0.1;
P.al1 = -0.95;
P.ga1 = k1/(U^2*L^3);
P.ga2 = k2/(U^2*L^3);
P.ze1 = 0.05;

% circuit (example)
P.lam   = 0.18;
P.kap_e = 0.2;
P.kap_c = 0.5;
P.sigma = 1.0;

sysP = [P.be1; P.be2; P.mu; P.al1; P.ga1; P.ze1; P.lam; ...
        P.kap_e; P.kap_c; P.sigma; P.ga2];

theta = sqrt(max(P.lam,0));

fprintf('================================================\n');
fprintf('Verify Energy/Consistency @ Omega=%.4f\n', Omega_test);
fprintf('Fw=%.4g | lam=%.4g kap_e=%.4g kap_c=%.4g sigma=%.4g | theta=%.4g\n', ...
    Fw, P.lam, P.kap_e, P.kap_c, P.sigma, theta);
fprintf('================================================\n');

%% ---------- 1) Initial guess from FRF ----------
fprintf('\n[1] Running FRF(sysP) for initial guess...\n');
x_frf = FRF(sysP);
if isempty(x_frf) || size(x_frf,1) < 16
    error('FRF output invalid.');
end
Om = x_frf(16,:).';
[~, idx0] = min(abs(Om - Omega_test));
y0_15 = x_frf(1:15, idx0);
fprintf('    Picked FRF point: Omega=%.6f (idx=%d)\n', Om(idx0), idx0);

%% ---------- 2) Solve fixed-Omega HBM by fsolve ----------
fprintf('\n[2] Solving fixed-Omega HBM by fsolve...\n');
fun = @(y15) nondim_temp2([y15(:); Omega_test], sysP);

opt = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxIterations',800, ...
    'MaxFunctionEvaluations',40000);

[y15_sol, fval, exitflag] = fsolve(fun, y0_15, opt);
fprintf('fsolve exitflag=%d, ||res||_inf=%.3e\n', exitflag, norm(fval,inf));

%% ---------- 3) Reconstruct time histories (HB 0,1,3) ----------
W = Omega_test;
T = 2*pi/W;
Nt = 6000;                    % DENSE sampling improves power balance accuracy
t  = linspace(0, T, Nt).';     % nondim time tau

x1c = y15_sol(1:5);
x2c = y15_sol(6:10);
qc  = y15_sol(11:15);

x1   = recon_u  (x1c, t, W);
x2   = recon_u  (x2c, t, W);
q    = recon_u  (qc , t, W);

x1p  = recon_up (x1c, t, W);
x2p  = recon_up (x2c, t, W);
qp   = recon_up (qc , t, W);        % q' (current, nondim)

x1pp = recon_upp(x1c, t, W);
qpp  = recon_upp(qc , t, W);

x12  = x1 - x2;
x12p = x1p - x2p;

Fexc = Fw*cos(W*t);

%% ---------- 4) HARD check A: Circuit residual in time domain ----------
% 【修正】：物理上反电动势是 +theta*(x1'-x2')，移项到等式左边做残差应为负号！
% R3(t) = kap_e*q'' + sigma*q' + kap_c*q - theta*(x1' - x2')
R3_t = P.kap_e*qpp + P.sigma*qp + P.kap_c*q - theta*x12p;

r3_rms_rel = rms(R3_t) / max(1e-14, rms(theta*x12p)) * 100;
fprintf('\n--- Circuit equation check (time domain) ---\n');
fprintf('RMS(R3) relative = %.6f %%   (R3 should be ~0)\n', r3_rms_rel);

%% ---------- 5) HARD check B: EM force consistency (HB-consistent cubic) ----------
% From your R1 structure:
%   x1'' + (be1+al1)*x12 + ga1*(x12^3) + theta*q' = Fexc
% => fem_LHS = Fexc - [x1'' + (be1+al1)*x12 + ga1*(x12^3)]
% => fem_RHS = theta*q'
f12 = (P.be1 + P.al1)*x12 + P.ga1*(x12.^3);

fem_LHS = Fexc - (x1pp + f12);
fem_RHS = theta*qp;

fem_rms_rel = rms(fem_LHS - fem_RHS) / max(1e-14, rms(fem_RHS)) * 100;
fprintf('\n--- EM force check (HB-consistent cubic) ---\n');
fprintf('RMS error = %.6f %%\n', fem_rms_rel);

%% ---------- 6) Cycle-averaged power balance ----------
% Input power
Pin = Fexc .* x1p;

% Mechanical dissipation (your model: 2*mu*ze1*x2' in eq => power = 2*mu*ze1*x2'^2)
Pmech = (2*P.mu*P.ze1) * (x2p.^2);

% Resistor dissipation (sigma*q' term => power = sigma*q'^2)
Pres  = P.sigma * (qp.^2);

Mean_Pin  = mean(Pin);
Mean_Pdiss = mean(Pmech + Pres);
Bal_err = abs(Mean_Pin - Mean_Pdiss) / max(1e-14, abs(Mean_Pin)) * 100;

fprintf('\n--- Cycle-averaged power balance ---\n');
fprintf('<Pin>    = %.8e\n', Mean_Pin);
fprintf('<Pdiss>  = %.8e   (= <Pmech> + <Pres>)\n', Mean_Pdiss);
fprintf('Balance err = %.3f %%\n', Bal_err);

%% ---------- 7) Plots ----------
figure('Color','w','Position',[180,80,980,860]);

subplot(3,1,1);
plot(t, x1, 'b', 'LineWidth',1.2); hold on;
plot(t, x2, 'r--', 'LineWidth',1.2);
grid on; xlim([0,T]);
ylabel('x');
legend('x_1','x_2');
title(sprintf('Time histories (\\Omega=%.3f, \\lambda=%.3f)', W, P.lam));

subplot(3,1,2);
plot(t, Pin, 'k', 'LineWidth',1.0); hold on;
plot(t, Pmech, 'r-.', 'LineWidth',1.0);
plot(t, Pres,  'g', 'LineWidth',1.2);
grid on; xlim([0,T]);
ylabel('Power');
legend('P_{in}','P_{mech}','P_{res}');
title(sprintf('Power terms (mean balance err = %.3f%%)', Bal_err));

subplot(3,1,3);
plot(t, R3_t, 'm', 'LineWidth',1.1);
grid on; xlim([0,T]);
xlabel('\tau'); ylabel('R_3(t)');
title(sprintf('Circuit residual R_3(t), RMS rel = %.4g%%', r3_rms_rel));

figure('Color','w','Position',[1180,160,860,420]);
plot(t, fem_LHS, 'k', 'LineWidth',1.1); hold on;
plot(t, fem_RHS, 'g--', 'LineWidth',1.1);
grid on; xlim([0,T]);
xlabel('\tau'); ylabel('F_{em}');
legend('F_{em,LHS} from mech eq','F_{em,RHS} = \\theta q''','Location','best');
title(sprintf('EM force consistency RMS rel = %.4g%%', fem_rms_rel));

%% -------- helper funcs (HB 0,1,3) --------
function u = recon_u(cfs, t, W)
    u = cfs(1) + cfs(2)*cos(W*t) + cfs(3)*sin(W*t) + cfs(4)*cos(3*W*t) + cfs(5)*sin(3*W*t);
end

function up = recon_up(cfs, t, W)
    up = cfs(2)*(-W*sin(W*t)) + cfs(3)*(W*cos(W*t)) + cfs(4)*(-3*W*sin(3*W*t)) + cfs(5)*(3*W*cos(3*W*t));
end

function upp = recon_upp(cfs, t, W)
    upp = cfs(2)*(-W^2*cos(W*t)) + cfs(3)*(-W^2*sin(W*t)) + cfs(4)*(-(3*W)^2*cos(3*W*t)) + cfs(5)*(-(3*W)^2*sin(3*W*t));
end