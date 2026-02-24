clc; clear; close all;

%% =========================================================
%  Strict H∞ optimization for resonance peak minimization
%  Outer variables: [sigma, kap_e, kap_c]
%  Inner loop: sweep frequency -> HBM(Newton) -> Floquet -> TF peak
%
%  Requires:
%     newton.m
%     nondim_temp2.m  (HB residual, uses global FixedOmega, Fw)
%
%  Notes:
%  - Force transmissibility is computed as:
%       TF(Ω) = |F_trans,1(Ω)| / Fw
%    where F_trans(t) = be2*x2 + ga2*x2^3 + 2*mu*ze*v2   (base transmitted force proxy)
%    and |·,1| is the 1st-harmonic amplitude via FFT.
%
%  - Floquet stability: compute max|mu| along HB orbit (Φ' = A(t)Φ)
%
%  You can change:
%    - frequency band [Om_min, Om_max]
%    - target sub-band (peak-band) if you only want to suppress one peak
%    - force level Fw_opt (optimization force level)
%
%% =========================================================

%% --------------------- Fixed physical parameters ---------------------
mu   = 0.2;
beta = 2.0;
K1   = 1.0;
K2   = 0.2;      % <- you used 0.2 here
U    = 2.0;
L    = 4/9;

v = 2.5;
alpha1 = v    - 2*K1*(1-L)/L;
alpha2 = beta - 2*K2*(1-L)/L;
gamma1 = K1/(U^2 * L^3);
gamma2 = K2/(U^2 * L^3);

P.be1 = 1.0;
P.al1 = alpha1 - P.be1;
P.be2 = alpha2;
P.ga1 = gamma1;
P.ga2 = gamma2;
P.mu  = mu;
P.ze1 = 0.05;       % used in x2 damping term in your ext_ode

P.lam = 0.18;       % theta^2
theta = sqrt(max(P.lam,0));

% Build sysP template (sigma/kap_e/kap_c will be overwritten by optimizer)
% sysP = [be1, be2, mu, al1, ga1, ze, lam, kap_e, kap_c, sigma, ga2]
sysP0 = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
         1.0,  0.2,  1.0,   P.ga2];

%% --------------------- Optimization settings -------------------------
% Force level used for frequency response & peak suppression
Fw_opt = 0.008;              % choose a representative force (can be your target region)

% Frequency sweep band (log sweep is common)
Om_min = 0.2;
Om_max = 6.0;
Nw     = 220;
Om_grid = logspace(log10(Om_min), log10(Om_max), Nw).';

% If you only want to suppress ONE peak: set a sub-band here, else use full band
use_peak_band = false;
Om_peak_center = 0.25;
peak_bw = 0.20;  % ±20%
idxBand = true(size(Om_grid));
if use_peak_band
    idxBand = (Om_grid >= Om_peak_center*(1-peak_bw)) & (Om_grid <= Om_peak_center*(1+peak_bw));
end

% Floquet settings
Nt_floquet = 600;             % RK4 steps per period for Phi
tol_stable = 1.002;

% Time grid for transmitted force FFT (AFT)
Nt_fft = 1024;                % points per period for FFT of F_trans
% (use power-of-two for speed)

% Outer search bounds for [sigma, kap_e, kap_c]
lb = [0.02, 0.02, 0.02];
ub = [3.00, 3.00, 3.00];

% Outer search: random + local refine
Nsamp = 800;       % random samples (increase to 1000+ if needed)
TopK  = 10;        % keep top-K for local refine

rng(1);

%% --------------------- Evaluate random candidates --------------------
X = lhsdesign(Nsamp,3);
X = lb + X.*(ub-lb);

J = inf(Nsamp,1);
cache = cell(Nsamp,1);

fprintf('Random screening: %d candidates...\n', Nsamp);

for i = 1:Nsamp
    sigma = X(i,1); kap_e = X(i,2); kap_c = X(i,3);

    [Ji, out] = objective_peak_Hinf_strict( ...
        sigma, kap_e, kap_c, sysP0, Om_grid, idxBand, Fw_opt, ...
        Nt_fft, Nt_floquet, tol_stable);

    J(i) = Ji;
    cache{i} = out;

    if mod(i,25)==0
        fprintf('  %4d/%4d  bestJ=%.4e\n', i, Nsamp, min(J));
    end
end

[Js, idx] = sort(J,'ascend');
idx = idx(isfinite(Js));
idx = idx(1:min(TopK,numel(idx)));

fprintf('\nTop-%d from random screening:\n', numel(idx));
fprintf(' rank |     J     |  sigma   kap_e   kap_c\n');
fprintf('------------------------------------------\n');
for r = 1:numel(idx)
    i = idx(r);
    fprintf('%4d | %8.3e | %6.3f  %6.3f  %6.3f\n', r, J(i), X(i,1), X(i,2), X(i,3));
end

%% --------------------- Local refinement (bounded fminsearch) ----------
% Use parameter transform to handle bounds: p = lb + (ub-lb).*sigmoid(z)
bestJ = inf; bestP = []; bestOut = [];

for r = 1:numel(idx)
    p0 = X(idx(r),:);

    z0 = inv_sigmoid((p0 - lb)./(ub-lb));
    funz = @(z) objective_z(z, lb, ub, sysP0, Om_grid, idxBand, Fw_opt, Nt_fft, Nt_floquet, tol_stable);

    opts = optimset('Display','iter','MaxIter',60,'TolX',1e-2,'TolFun',1e-3);
    [zopt, Jopt] = fminsearch(funz, z0, opts);

    popt = lb + (ub-lb).*sigmoid(zopt);

    [Jtrue, out] = objective_peak_Hinf_strict( ...
        popt(1), popt(2), popt(3), sysP0, Om_grid, idxBand, Fw_opt, ...
        Nt_fft, Nt_floquet, tol_stable);

    fprintf('\nLocal refine #%d done: J=%.4e, p=[%.4f %.4f %.4f]\n', r, Jtrue, popt);

    if Jtrue < bestJ
        bestJ = Jtrue;
        bestP = popt;
        bestOut = out;
    end
end

fprintf('\n================== BEST RESULT ==================\n');
fprintf('Best J(H∞ peak + stability) = %.6e\n', bestJ);
fprintf('Best [sigma, kap_e, kap_c]  = [%.6f, %.6f, %.6f]\n', bestP(1), bestP(2), bestP(3));
fprintf('=================================================\n');

%% --------------------- Plot best candidate response -------------------
figure('Color','w','Position',[120 120 1100 420]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; box on; grid on;
loglog(Om_grid, bestOut.TF, 'LineWidth',1.7);
xlabel('\Omega'); ylabel('T_F(\Omega)=|F_{t,1}|/F_w');
title('Force transmissibility (best)');
if use_peak_band
    xline(Om_peak_center*(1-peak_bw),'k--');
    xline(Om_peak_center*(1+peak_bw),'k--');
end

nexttile; hold on; box on; grid on;
semilogx(Om_grid, bestOut.maxMu, 'LineWidth',1.7);
yline(tol_stable,'--','LineWidth',1.4);
xlabel('\Omega'); ylabel('max|\mu|');
title('Floquet stability (best)');
ylim([0, max(1.5, 1.05*max(bestOut.maxMu))]);

%% =========================================================
%  Nested functions
%% =========================================================

function J = objective_z(z, lb, ub, sysP0, Om_grid, idxBand, Fw, Nt_fft, Nt_floquet, tol_stable)
    p = lb + (ub-lb).*sigmoid(z);
    [J, ~] = objective_peak_Hinf_strict(p(1),p(2),p(3), sysP0, Om_grid, idxBand, Fw, Nt_fft, Nt_floquet, tol_stable);
end

function s = sigmoid(z)
    s = 1./(1+exp(-z));
end

function z = inv_sigmoid(s)
    s = min(max(s,1e-6), 1-1e-6);
    z = log(s./(1-s));
end

% end main script


%% =========================================================
%  Objective: strict peak minimization + Floquet stability
%% =========================================================
function [J, out] = objective_peak_Hinf_strict( ...
    sigma, kap_e, kap_c, sysP0, Om_grid, idxBand, Fw, Nt_fft, Nt_floquet, tol_stable)

    % Build sysP
    sysP = sysP0;
    sysP(8)  = kap_e;
    sysP(9)  = kap_c;
    sysP(10) = sigma;

    % Hard reject: avoid near singular operator (den close to 0)
    den = kap_e*Om_grid.^2 - 1i*sigma*Om_grid - kap_c;
    if min(abs(den)) < 1e-2
        J = 1e6 + 1e3*(1e-2/min(abs(den)));
        out = empty_out(Om_grid);
        return;
    end

    % Inner sweep using HBM/Newton:
    % Use continuation by carrying previous HB solution as initial guess.
    global FixedOmega Fw_global
    Fw_global = Fw;

    Nw = numel(Om_grid);
    TF = nan(Nw,1);
    maxMu = nan(Nw,1);
    ok = false(Nw,1);

    % initial guess: zeros
    y_guess = [zeros(15,1); Fw];

    fail_count = 0;

    for j = 1:Nw
        Om = Om_grid(j);
        FixedOmega = Om;

        % Solve HB residual
        try
            y_sol = newton('nondim_temp2', y_guess, sysP);
        catch
            fail_count = fail_count + 1;
            if fail_count > 8
                break;
            else
                continue;
            end
        end

        x_coeff = y_sol(1:15);
        y_guess = [x_coeff; Fw];  % continuation initial guess

        ok(j) = true;

        % Compute transmitted force TF at this Omega
        TF(j) = compute_force_transmissibility_fromHB(x_coeff, sysP, Om, Fw, Nt_fft);

        % Floquet stability at this Omega
        maxMu(j) = compute_floquet_maxmu_fromHB(x_coeff, sysP, Om, Fw, Nt_floquet);

        % Early penalty if clearly unstable
        if maxMu(j) > 1.10
            % still continue a bit is ok, but we can break to save time
            % break;
        end
    end

    % If too many failures, penalize
    if nnz(ok) < 0.7*Nw
        J = 1e6 + 1e3*(1 - nnz(ok)/Nw);
        out = pack_out(Om_grid, TF, maxMu);
        return;
    end

    % Only evaluate band
    TFb = TF(idxBand & ok);
    mub = maxMu(idxBand & ok);

    if isempty(TFb)
        J = 1e6;
        out = pack_out(Om_grid, TF, maxMu);
        return;
    end

    % Peak objective: H∞
    peakTF = max(TFb);

    % Stability penalty: if any point violates tol, add strong penalty
    vio = max(0, mub - tol_stable);
    pen_stab = 0;
    if any(vio > 0)
        pen_stab = 5e2 * max(vio)^2 + 2e2 * mean(vio > 0);
    end

    % Small smoothness regularizer (optional)
    pen_nan = 0;
    if any(~isfinite(TFb)) || any(~isfinite(mub))
        pen_nan = 1e3;
    end

    J = peakTF + pen_stab + pen_nan;

    out = pack_out(Om_grid, TF, maxMu);
end

function out = empty_out(Om)
    out = struct();
    out.TF = nan(size(Om));
    out.maxMu = nan(size(Om));
end

function out = pack_out(Om, TF, maxMu)
    out = struct();
    out.TF = TF;
    out.maxMu = maxMu;
end

%% =========================================================
%  Force transmissibility from HB coefficients (AFT + FFT)
%  F_trans(t) = be2*x2 + ga2*x2^3 + 2*mu*ze*v2
%  TF = |F_trans,1| / Fw
%% =========================================================
function TF = compute_force_transmissibility_fromHB(x_coeff, sysP, Omega, Fw, Nt)
    % unpack
    x1c = x_coeff(1:5);
    x2c = x_coeff(6:10);
    qc  = x_coeff(11:15); %#ok<NASGU>

    be2 = sysP(2);
    mu  = sysP(3);
    ze  = sysP(6);
    ga2 = sysP(11);

    % time grid over one period
    T = 2*pi/Omega;
    t = linspace(0, T, Nt+1); t(end) = [];
    % reconstruct x2, v2
    [x2, v2] = reconstruct_onevar(t, x2c, Omega);

    % transmitted force proxy
    Ftr = be2*x2 + ga2*(x2.^3) + 2*mu*ze*v2;

    % FFT first harmonic amplitude
    Y = fft(Ftr)/Nt;
    % index 2 corresponds to 1-cycle harmonic when t spans exactly one period
    A1 = 2*abs(Y(2)); % single-sided amplitude

    TF = A1 / max(Fw, 1e-12);
end

%% =========================================================
%  Floquet max|mu| from HB coefficients (HB-orbit based)
%  Φ' = A(t)Φ , Φ(0)=I, max|eig(Φ(T))|
%% =========================================================
function maxMu = compute_floquet_maxmu_fromHB(x_coeff, sysP, Omega, Fw, Nt)
    %#ok<NASGU> Fw not used in A(t) because forcing does not enter Jacobian

    x1c = x_coeff(1:5);
    x2c = x_coeff(6:10);
    qc  = x_coeff(11:15);

    T = 2*pi/Omega;
    dt = T/Nt;

    Phi = eye(6);
    t = 0;

    for i = 1:Nt
        k1 = dPhi_full(t,          Phi, sysP, Omega, x1c, x2c, qc);
        k2 = dPhi_full(t + 0.5*dt, Phi + 0.5*dt*k1, sysP, Omega, x1c, x2c, qc);
        k3 = dPhi_full(t + 0.5*dt, Phi + 0.5*dt*k2, sysP, Omega, x1c, x2c, qc);
        k4 = dPhi_full(t + dt,     Phi + dt*k3,     sysP, Omega, x1c, x2c, qc);
        Phi = Phi + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end

    ev = eig(Phi);
    maxMu = max(abs(ev));
end

function dPhi = dPhi_full(t, Phi, sysP, Omega, x1c, x2c, qc)
    % reconstruct state along HB orbit
    y = HB_reconstruct_state(t, x1c, x2c, qc, Omega);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));

    x1=y(1); v1=y(2); x2=y(3); v2=y(4); q=y(5); qd=y(6); %#ok<NASGU>
    dx = x1-x2;

    df12_ddx = (be1+al1) + 3*ga1*dx^2;
    df2g_dx2 = be2 + 3*ga2*x2^2;
    df2g_dv2 = 2*mu*ze;

    A = zeros(6);
    A(1,2) = 1;
    A(2,1) = -df12_ddx;      A(2,3) = +df12_ddx;                      A(2,6) = theta;
    A(3,4) = 1;
    A(4,1) = df12_ddx/mu;    A(4,3) = (-df12_ddx - df2g_dx2)/mu;      A(4,4) = -df2g_dv2/mu;   A(4,6) = -theta/mu;
    A(5,6) = 1;

    A(6,2) = -theta/kap_e;
    A(6,4) = +theta/kap_e;
    A(6,5) = -kap_c/kap_e;
    A(6,6) = -sigma/kap_e;

    dPhi = A*Phi;
end

function y = HB_reconstruct_state(t, x1c, x2c, qc, Omega)
    [x1, v1] = reconstruct_onevar(t, x1c, Omega);
    [x2, v2] = reconstruct_onevar(t, x2c, Omega);
    [q,  qd] = reconstruct_onevar(t, qc,  Omega);
    y = [x1(:).'; v1(:).'; x2(:).'; v2(:).'; q(:).'; qd(:).'];
    y = y(:); % 6x1 if t is scalar, or 6*Nt stacked if vector (we use scalar here)
end

function [x, v] = reconstruct_onevar(t, c, Omega)
    % c = [a0 a1c a1s a3c a3s]
    w = Omega;
    ct  = cos(w*t);  st  = sin(w*t);
    c3t = cos(3*w*t); s3t = sin(3*w*t);

    x = c(1) + c(2)*ct + c(3)*st + c(4)*c3t + c(5)*s3t;
    v = (-w*c(2))*st + (w*c(3))*ct + (-3*w*c(4))*s3t + (3*w*c(5))*c3t;
end