%% Run_BoA_Space_FIXED.m
% Basins of Attraction (BoA) for 6D electromechanical 2DOF system
% Slice: (x2(0), v2(0)) grid, other states fixed at a reference phase point.
% Classification: nearest attractor (LOW vs HIGH) using strict Poincare anchors (kT).
%
% Key FIX:
%   - anchors are extracted by FIXED-STEP RK4 and sampled exactly at kT (no ode45 phase drift)
%   - add diverged/unknown class
%   - sanity checks for anchor separation

clear; clc; close all;

%% =========================
% 1) System params
% =========================
% sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
P.be1 = 1.0;  P.be2 = 0.1;  P.mu  = 0.2;  P.al1 = -0.95;
P.ga1 = 1.5;  P.ga2 = 1.5;  P.ze1 = 0.05;
P.lam = 0.18; P.kap_e = 1.0; P.kap_c = 0.2; P.sigma = 0.5;

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, ...
        P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

Omega_test = 0.25;
Fw_test    = 0.00843978;   % 用你验证双稳态的 MID 点

%% =========================
% 2) Provide TWO ICs that land on LOW/HIGH (IMPORTANT!)
% =========================
% 这里必须换成你“专门对 MID 做双初值积分对照”时用的两个初值
% Format: [x1; v1; x2; v2; q; qd]
y0_low  = [0;0; -0.02; 0.00; 0; 0];   % <-- 替换成你的 LOW IC
y0_high = [0;0;  0.20; 0.00; 0; 0];   % <-- 替换成你的 HIGH IC

%% =========================
% 3) BoA slice settings
% =========================
N_grid  = 250;
x2_span = linspace(-0.25, 0.25, N_grid);
v2_span = linspace(-0.15, 0.15, N_grid);
[X2_0, V2_0] = meshgrid(x2_span, v2_span);

slice_ref = "low";   % "low" or "high" : fix other 4 states at which attractor's phase point

% weights for distance in 6D (建议先用全1；如果电变量量级差太大再调)
w = [1, 1, 1, 1, 1, 1];
W = diag(w);

%% =========================
% 4) Integration controls (RK4)
% =========================
N_per_cycle   = 100;    % 建议 80~160；你非线性强可以 120/160
N_settle_cyc  = 250;    % anchor settle cycles
N_keep        = 15;     % keep last N Poincare points for averaging anchor and classification
N_total_cyc   = 180;    % BoA total cycles
N_drop_cyc    = 140;    % drop transient cycles

% diverge detection
blow_limit = 50;        % 状态超过这个认为发散（无量纲下足够大）
nan_limit  = 1e6;

%% =========================
% 5) Get strict Poincare anchors by RK4 (kT samples)
% =========================
fprintf('Step 1/3: Extracting LOW/HIGH anchors using RK4 @ exact kT...\n');

[anchor_low,  rms_low]  = get_anchor_rk4(sysP, Omega_test, Fw_test, y0_low,  N_settle_cyc, N_per_cycle, N_keep, blow_limit);
[anchor_high, rms_high] = get_anchor_rk4(sysP, Omega_test, Fw_test, y0_high, N_settle_cyc, N_per_cycle, N_keep, blow_limit);

fprintf('  LOW  anchor:  [%.4g %.4g %.4g %.4g %.4g %.4g] | rms(x2)=%.4g\n', anchor_low,  rms_low);
fprintf('  HIGH anchor:  [%.4g %.4g %.4g %.4g %.4g %.4g] | rms(x2)=%.4g\n', anchor_high, rms_high);

if norm(anchor_low - anchor_high) < 1e-2 && abs(rms_low - rms_high) < 1e-2
    warning(['LOW/HIGH anchors look too close. ' ...
             'Likely your two ICs landed on the SAME attractor (or F is not bistable).']);
end

% choose slicing phase point
if slice_ref == "low"
    y_fix = anchor_low(:);
else
    y_fix = anchor_high(:);
end

%% =========================
% 6) Build initial matrix for BoA slice
% =========================
fprintf('Step 2/3: Building slice ICs...\n');
M_points = N_grid * N_grid;
Y = zeros(6, M_points);

% fixed dims at chosen phase point
Y(1,:) = y_fix(1);  % x1
Y(2,:) = y_fix(2);  % v1
Y(5,:) = y_fix(5);  % q
Y(6,:) = y_fix(6);  % qd

% scanned dims
Y(3,:) = X2_0(:).';
Y(4,:) = V2_0(:).';

%% =========================
% 7) Vectorized RK4 integrate and sample last N_keep Poincare points
% =========================
fprintf('Step 3/3: Vectorized RK4 BoA integration (%d points)...\n', M_points);

T  = 2*pi/Omega_test;
dt = T / N_per_cycle;

keep_buf = zeros(6, M_points, N_keep);
keep_idx = 0;

diverged = false(1, M_points);

tic;
for cyc = 1:N_total_cyc
    for st = 1:N_per_cycle
        t = (cyc-1)*T + (st-1)*dt;

        k1 = sys_ode_vec(t,          Y,               sysP, Omega_test, Fw_test);
        k2 = sys_ode_vec(t + 0.5*dt, Y + 0.5*dt.*k1,  sysP, Omega_test, Fw_test);
        k3 = sys_ode_vec(t + 0.5*dt, Y + 0.5*dt.*k2,  sysP, Omega_test, Fw_test);
        k4 = sys_ode_vec(t + dt,     Y + dt.*k3,      sysP, Omega_test, Fw_test);
        Y  = Y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

        % divergence check (vectorized)
        bad = any(~isfinite(Y),1) | any(abs(Y) > blow_limit,1) | any(abs(Y) > nan_limit,1);
        diverged = diverged | bad;

        % for diverged points, clamp to avoid NaN spreading
        if any(bad)
            Y(:,bad) = 0;
        end
    end

    % Poincare sample at exact cycle boundary
    if cyc > (N_total_cyc - N_keep)
        keep_idx = keep_idx + 1;
        keep_buf(:,:,keep_idx) = Y;
    end

    if mod(cyc, 30)==0
        fprintf('  cycles %d / %d done\n', cyc, N_total_cyc);
    end
end
toc;

Y_end = mean(keep_buf, 3); % averaged Poincare state

%% =========================
% 8) Classification (LOW/HIGH/Diverged)
% =========================
d_low  = distW(Y_end, anchor_low(:),  W);
d_high = distW(Y_end, anchor_high(:), W);

% BoA class:
%   0 = LOW, 1 = HIGH, 2 = Diverged/unknown
BoA = zeros(1, M_points);
BoA(d_high < d_low) = 1;
BoA(diverged) = 2;

ratio_low  = 100*mean(BoA==0);
ratio_high = 100*mean(BoA==1);
ratio_div  = 100*mean(BoA==2);
fprintf('\nBoA ratio: LOW=%.2f%%, HIGH=%.2f%%, DIV=%.2f%% (slice_ref=%s)\n', ...
    ratio_low, ratio_high, ratio_div, slice_ref);

BoA_Map = reshape(BoA, N_grid, N_grid);

%% =========================
% 9) Plot (SCI-style)
% =========================
figure('Color','w','Position',[120 80 860 680]);
ax = gca; hold on; box on;

imagesc(x2_span, v2_span, BoA_Map);
set(ax,'YDir','normal');

% colormap: 0 LOW blue, 1 HIGH red, 2 DIV gray
cmap = [0, 0.4470, 0.7410;
        0.8500, 0.3250, 0.0980;
        0.65, 0.65, 0.65];
colormap(cmap);
caxis([0 2]);

% boundary between low/high: contour of (d_high - d_low)=0 (only meaningful where not diverged)
Ddiff = reshape(d_high - d_low, N_grid, N_grid);
contour(x2_span, v2_span, Ddiff, [0 0], 'k', 'LineWidth', 1.1);

% mark slicing fixed point projection (x2,v2)
plot(y_fix(3), y_fix(4), 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor','y', 'MarkerEdgeColor','k', 'LineWidth', 1.3);

xlabel('Initial displacement  x_2(0)','FontName','Times New Roman','FontSize',14);
ylabel('Initial velocity      v_2(0)','FontName','Times New Roman','FontSize',14);
title(sprintf('Basins of attraction (nearest-attractor)  \\Omega=%.3g, F=%.5g',Omega_test,Fw_test), ...
      'FontName','Times New Roman','FontSize',14);

h1 = patch(NaN,NaN,cmap(1,:),'EdgeColor','k');
h2 = patch(NaN,NaN,cmap(2,:),'EdgeColor','k');
h3 = patch(NaN,NaN,cmap(3,:),'EdgeColor','k');
legend([h1 h2 h3], {'LOW attractor (safe)','HIGH attractor (dangerous)','Diverged/unknown'}, ...
    'Location','northoutside','Orientation','horizontal', ...
    'FontName','Times New Roman','FontSize',12);

set(gca,'FontName','Times New Roman','FontSize',12);
axis square; grid on;

%% =========================
% ======== local functions ========
% =========================
function [anchor, rms_x2] = get_anchor_rk4(sysP, Omega, Fw, y0, N_settle_cyc, N_per_cycle, N_keep, blow_limit)
    % settle by RK4 and sample exact Poincare points at kT
    T  = 2*pi/Omega;
    dt = T/N_per_cycle;

    y = y0(:);
    keep = zeros(6, N_keep);
    keep_idx = 0;

    for cyc = 1:N_settle_cyc
        for st = 1:N_per_cycle
            t = (cyc-1)*T + (st-1)*dt;

            k1 = sys_ode_scalar(t,          y,               sysP, Omega, Fw);
            k2 = sys_ode_scalar(t + 0.5*dt, y + 0.5*dt*k1,    sysP, Omega, Fw);
            k3 = sys_ode_scalar(t + 0.5*dt, y + 0.5*dt*k2,    sysP, Omega, Fw);
            k4 = sys_ode_scalar(t + dt,     y + dt*k3,        sysP, Omega, Fw);
            y  = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

            if any(~isfinite(y)) || any(abs(y) > blow_limit)
                warning('Anchor settling diverged. Check IC / F / dt.');
                anchor = nan(1,6);
                rms_x2 = nan;
                return;
            end
        end

        % record last N_keep cycle-boundary points
        if cyc > (N_settle_cyc - N_keep)
            keep_idx = keep_idx + 1;
            keep(:,keep_idx) = y;
        end
    end

    % anchor = mean of last N_keep Poincare points (period-1 orbit -> almost identical)
    anchor = mean(keep,2).';
    x2_series = keep(3,:);
    rms_x2 = sqrt(mean(x2_series.^2));
end

function d = distW(Y, a, W)
    % weighted Euclidean distance column-wise
    % Y: 6xM, a: 6x1
    Z = Y - a;
    d = sqrt(sum((W*Z).^2, 1)).';
end

function dy = sys_ode_scalar(t, y, sysP, Omega, Fw)
    Y = y(:);
    dY = sys_ode_vec(t, Y, sysP, Omega, Fw);
    dy = dY(:);
end

function dY = sys_ode_vec(t, Y, sysP, Omega, Fw)
    % Vectorized ODE: Y can be 6xM or 6x1
    x1=Y(1,:); v1=Y(2,:); x2=Y(3,:); v2=Y(4,:); q=Y(5,:); qd=Y(6,:);

    be1=sysP(1); be2=sysP(2); mu=sysP(3);
    al1=sysP(4); ga1=sysP(5); ze2=sysP(6);
    lam=sysP(7); kap_e=sysP(8); kap_c=sysP(9); sigma=sysP(10); ga2=sysP(11);

    theta = sqrt(max(lam,0));

    dx = x1 - x2;
    dv = v1 - v2;

    f12 = (be1+al1).*dx + ga1.*(dx.^3);
    f2g = be2.*x2 + ga2.*(x2.^3) + 2*mu*ze2.*v2;

    % NOTE: your model uses +theta*qd in x1 eq and -theta*qd in x2 eq
    x1dd = -f12 + theta.*qd + Fw*cos(Omega*t);
    x2dd = (f12 - f2g - theta.*qd) / mu;

    if kap_e == 0
        s = max(abs(sigma), 1e-12);
        qd_new = (-kap_c.*q - theta.*dv) / s;
        qdd = (qd_new - qd) .* 50;
    else
        qdd = (-sigma.*qd - kap_c.*q - theta.*dv) / kap_e;
    end

    dY = [v1; x1dd; v2; x2dd; qd; qdd];
end
