%% Verify_Floquet_Stability_Current.m
% Floquet stability of a periodic solution (HBM fixed-Omega result)
% For your current model:
%   sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
% HBM vector y15:
%   1:5   x1  [x10, a11, b11, a13, b13]
%   6:10  x2  [x20, a21, b21, a23, b23]
%   11:15 q   [q0,  aq1, bq1, aq3, bq3]
%
% State for ODE/Floquet:
%   y = [x1; v1; x2; v2; q; qd]
%
% Stability criterion:
%   all |eig(Phi(T))| < 1  => stable periodic orbit
%
% NOTE:
%   This checks stability of the time-domain periodic orbit corresponding to HB solution.

clear; clc; close all;

%% -----------------------------
% 0) set case
% -----------------------------
global Fw
Fw = 0.005;

P.be1   = 1.0;
P.be2   = 0.1;
P.mu    = 0.2;
P.al1   = -0.95;
% Example nonlinear coeffs (replace by your k1,k2,L,U mapping if you want)
L = 4/9; U = 2; k1 = 1; k2 = 0.8;
P.ga1   = k1/(U^2*L^3);
P.ga2   = k2/(U^2*L^3);
P.ze1   = 0.05;

% circuit
P.lam   = 0.18;
P.kap_e = 0.2;
P.kap_c = 0.5;
P.sigma = 1.0;

sysP = [P.be1; P.be2; P.mu; P.al1; P.ga1; P.ze1; P.lam; ...
        P.kap_e; P.kap_c; P.sigma; P.ga2];

Omega = 0.80;

fprintf('====================================================\n');
fprintf('Floquet stability check (HBM->ODE variational)\n');
fprintf('Omega=%.4f, Fw=%.4f | lam=%.4f kap_e=%.4f kap_c=%.4f sigma=%.4f\n', ...
    Omega, Fw, P.lam, P.kap_e, P.kap_c, P.sigma);
fprintf('====================================================\n');

%% -----------------------------
% 1) Solve fixed-Omega HBM (fsolve) to get y_sol15
%    You can feed your own y_init15 if you already have it.
% -----------------------------
% Option A: use FRF point as initial guess (recommended if FRF works)
y_init15 = zeros(15,1);
try
    x_frf = FRF(sysP);
    if ~isempty(x_frf) && size(x_frf,1) >= 16
        Om = x_frf(16,:).';
        [~, idx] = min(abs(Om - Omega));
        y_init15 = x_frf(1:15, idx);
        fprintf('[Init] from FRF at Omega≈%.6f\n', Om(idx));
    else
        fprintf('[Init] FRF empty -> use small random init\n');
        y_init15 = 1e-3*randn(15,1);
    end
catch
    fprintf('[Init] FRF failed -> use small random init\n');
    y_init15 = 1e-3*randn(15,1);
end

solve_fun = @(y15) nondim_temp2([y15(:); Omega], sysP);

opt = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'MaxIterations',800, ...
    'MaxFunctionEvaluations',40000);

fprintf('\n[1] Solving fixed-Omega HBM by fsolve...\n');
[y_sol15, fval, exitflag] = fsolve(solve_fun, y_init15, opt);
fprintf('fsolve exitflag=%d, ||res||_inf=%.3e\n', exitflag, norm(fval,inf));

if exitflag <= 0
    error('HBM fixed-Omega solve failed. Try different init or smaller Fw.');
end

%% -----------------------------
% 2) Build initial state at t0=0 from HB coefficients
% -----------------------------
x1c = y_sol15(1:5);
x2c = y_sol15(6:10);
qc  = y_sol15(11:15);

t0 = 0;

x1  = recon_u (x1c, t0, Omega);
v1  = recon_up(x1c, t0, Omega);
x2  = recon_u (x2c, t0, Omega);
v2  = recon_up(x2c, t0, Omega);
q   = recon_u (qc,  t0, Omega);
qd  = recon_up(qc,  t0, Omega);

y0 = [x1; v1; x2; v2; q; qd];

%% -----------------------------
% 3) Integrate one period with variational equations
% -----------------------------
T = 2*pi/Omega;

% initial Phi = I
Phi0 = eye(6);
Y0_ext = [y0; Phi0(:)];

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

fprintf('\n[2] Integrating state+variational over one period...\n');
[tt, YY] = ode45(@(t,Y) ext_ode(t, Y, sysP, Omega, Fw), [0, T], Y0_ext, opts);

Y_end = YY(end,:).';
PhiT  = reshape(Y_end(7:end), 6,6);

mu = eig(PhiT);

% remove the trivial one near 1 (phase invariance) just for reporting
abs_mu = sort(abs(mu),'descend');
rho = max(abs(mu));

fprintf('\n================== Floquet Result ==================\n');
fprintf('Max |mu| = %.6f\n', rho);
fprintf('Eigenvalues |mu| (desc):\n');
disp(abs_mu.');

if rho < 1
    fprintf('[STABLE] periodic orbit is Floquet-stable at this Omega.\n');
else
    fprintf('[UNSTABLE] periodic orbit is Floquet-unstable at this Omega.\n');
end
fprintf('====================================================\n');

%% -----------------------------
% 4) Plot multipliers
% -----------------------------
figure('Color','w'); 
plot(real(mu), imag(mu), 'o','LineWidth',1.5); grid on; axis equal;
hold on;
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), '--');
xlabel('Re(\mu)'); ylabel('Im(\mu)');
title(sprintf('Floquet multipliers @ \\Omega=%.3f (max|\\mu|=%.3f)', Omega, rho));

%% ===== helper: extended ODE (state + variational) =====
function dY = ext_ode(t, Y, sysP, Omega, Fw)
    y   = Y(1:6);
    Phi = reshape(Y(7:end), 6,6);

    % state derivative
    dy  = odesys_nd(t, y, sysP, Omega, Fw);

    % Jacobian A(t) = df/dy
    A   = jac_nd(t, y, sysP, Omega, Fw);

    dPhi = A * Phi;

    dY = [dy; dPhi(:)];
end

%% ===== your current nondim ODE (consistent with your nondim_temp2 sign) =====
function dydt = odesys_nd(t, y, sysP, Omega, Fw)
    x1=y(1); v1=y(2);
    x2=y(3); v2=y(4);
    q =y(5); qd=y(6);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));

    dx  = (x1-x2);
    dv  = (v1-v2);

    % nonlinear coupling force (between masses)
    f12 = (be1+al1)*dx + ga1*dx^3;

    % lower-to-ground force (include damping on x2)
    f2g = be2*x2 + ga2*x2^3 + 2*mu*ze2*v2;

    % mechanical equations (match your nondim_temp2: +theta*qd in R1 and -theta*qd in R2)
    x1dd = -f12 + theta*qd + Fw*cos(Omega*t);
    x2dd = ( f12 - f2g - theta*qd )/mu;

    % circuit equation:
    % kap_e*qdd + sigma*qd + kap_c*q = - theta*(v1 - v2)
    % (THIS sign is critical for energy dissipation)
    if kap_e == 0
        tiny = 1e-12;
        qd_new = (-kap_c*q - theta*dv)/max(abs(sigma), tiny);
        qdd = (qd_new - qd)*50;
    else
        qdd = (-sigma*qd - kap_c*q - theta*dv)/kap_e;
    end

    dydt = [v1; x1dd; v2; x2dd; qd; qdd];
end

%% ===== Jacobian of ODE for variational equation =====
function A = jac_nd(t, y, sysP, Omega, Fw) %#ok<INUSD>
    x1=y(1); v1=y(2);
    x2=y(3); v2=y(4);
    q =y(5); qd=y(6); %#ok<NASGU>

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));

    dx = x1-x2;

    % partials
    df12_ddx = (be1+al1) + 3*ga1*dx^2;
    df2g_dx2 = be2 + 3*ga2*x2^2;
    df2g_dv2 = 2*mu*ze2;

    % state: [x1 v1 x2 v2 q qd]
    A = zeros(6);

    % x1' = v1
    A(1,2)=1;

    % v1' = x1dd = -f12 + theta*qd + Fw*cos(...)
    % d/dx1: -df12/ddx * d(dx)/dx1 = -df12_ddx
    % d/dx2: -df12/ddx * d(dx)/dx2 = +df12_ddx
    A(2,1) = -df12_ddx;
    A(2,3) = +df12_ddx;
    A(2,6) = +theta;

    % x2' = v2
    A(3,4)=1;

    % v2' = x2dd = ( f12 - f2g - theta*qd )/mu
    % d/dx1: +(df12_ddx)/mu
    % d/dx2: (-df12_ddx - df2g_dx2)/mu
    % d/dv2: (-df2g_dv2)/mu
    % d/dqd: (-theta)/mu
    A(4,1) = +df12_ddx/mu;
    A(4,3) = (-df12_ddx - df2g_dx2)/mu;
    A(4,4) = (-df2g_dv2)/mu;
    A(4,6) = (-theta)/mu;

    % q' = qd
    A(5,6)=1;

    % qd' = qdd
    % kap_e*qdd + sigma*qd + kap_c*q = -theta*(v1-v2)
    % => qdd = (-sigma*qd - kap_c*q - theta*(v1-v2))/kap_e
    if kap_e ~= 0
        A(6,2) = -theta/kap_e;      % v1
        A(6,4) = +theta/kap_e;      % v2
        A(6,5) = -kap_c/kap_e;      % q
        A(6,6) = -sigma/kap_e;      % qd
    else
        % for kap_e=0 relaxation form (approx)
        % qdd ≈ 50*(qd_new-qd), qd_new = (-kap_c*q - theta*(v1-v2))/sigma
        % => qdd = 50*( (-kap_c/sigma)*q + (-theta/sigma)*v1 + (theta/sigma)*v2 - qd )
        tiny=1e-12;
        s = max(abs(sigma),tiny);
        A(6,2)= 50*(-theta/s);
        A(6,4)= 50*(+theta/s);
        A(6,5)= 50*(-kap_c/s);
        A(6,6)= -50;
    end
end

%% ===== HB reconstruction helpers =====
function u = recon_u(cfs, t, W)
    u = cfs(1) + cfs(2)*cos(W*t) + cfs(3)*sin(W*t) + cfs(4)*cos(3*W*t) + cfs(5)*sin(3*W*t);
end
function up = recon_up(cfs, t, W)
    up = cfs(2)*(-W*sin(W*t)) + cfs(3)*(W*cos(W*t)) + cfs(4)*(-3*W*sin(3*W*t)) + cfs(5)*(3*W*cos(3*W*t));
end
