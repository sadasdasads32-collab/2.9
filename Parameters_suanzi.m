clc; clear; close all;
%三个算子的分布
%% ---- Parameters ----
theta = 0.18;
kap_e = 2.00;
kap_c = 0.15;
sigma = 0.10;

Om_min = 0.1;
Om_max = 10;
Npts   = 800;
Om = logspace(log10(Om_min), log10(Om_max), Npts).';

%% ---- Core operator ----
den = kap_e*Om.^2 - 1i*sigma*Om - kap_c;
K   = (theta^2).*Om.^2 ./ den;

Kr = real(K);
Ki = imag(K);

Ceq = Ki ./ max(Om,1e-12);

Om_e = sqrt(kap_c/kap_e);
fprintf('Omega_e = sqrt(kap_c/kap_e) = %.6f\n', Om_e);
fprintf('K_HF = theta^2/kap_e = %.6f\n', theta^2/kap_e);
fprintf('M_LF = theta^2/kap_c = %.6f\n', theta^2/kap_c);

%% ==============================================================
%  Split A: HF-anchored stiffness (constant Keq)
% ==============================================================
Keq_A = (theta^2)/kap_e * ones(size(Om));        % constant
Meq_A = (Keq_A - Kr) ./ max(Om.^2,1e-12);        % residual inertia

%% ==============================================================
%  Split B (RECOMMENDED): LF-anchored inertia (constant Meq)
%  Meq fixed to LF limit; Keq becomes frequency dependent
% ==============================================================
Meq_B = (theta^2)/kap_c * ones(size(Om));        % constant inertia
Keq_B = Kr + Om.^2 .* Meq_B;                     % stiffness compensation

%% ---- Plot 1: Kr & Ki ----
figure('Color','w','Position',[100 100 900 420]);
semilogx(Om, Kr, 'LineWidth',1.6); hold on;
semilogx(Om, Ki, 'LineWidth',1.6);
grid on; box on;
xline(Om_e,'k--','\Omega_e','LabelOrientation','horizontal');
xlabel('\Omega (log scale)'); ylabel('\mathcal{K}_r(\Omega), \mathcal{K}_i(\Omega)');
title('Real/Imag parts of \mathcal{K}(\Omega)');
legend('\mathcal{K}_r','\mathcal{K}_i','Location','best');

%% ---- Plot 2: Ceq ----
figure('Color','w','Position',[130 130 900 360]);
semilogx(Om, Ceq, 'LineWidth',1.8);
grid on; box on;
xline(Om_e,'k--','\Omega_e','LabelOrientation','horizontal');
xlabel('\Omega (log scale)'); ylabel('C_{eq}(\Omega)=Im(\mathcal{K})/\Omega');
title('Equivalent damping (unique definition)');

%% ---- Plot 3: Compare splits ----
figure('Color','w','Position',[160 160 1120 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

% Inertia
nexttile;
semilogx(Om, Meq_A, 'LineWidth',1.5); hold on;
semilogx(Om, Meq_B, 'LineWidth',1.5);
grid on; box on;
xline(Om_e,'k--');
xlabel('\Omega'); ylabel('M_{eq}(\Omega)');
title('Equivalent inertia');
legend('Split A: HF-anchored','Split B: LF-anchored','Location','best');

% Stiffness
nexttile;
semilogx(Om, Keq_A, 'LineWidth',1.5); hold on;
semilogx(Om, Keq_B, 'LineWidth',1.5);
grid on; box on;
xline(Om_e,'k--');
xlabel('\Omega'); ylabel('K_{eq}(\Omega)');
title('Equivalent stiffness');
legend('Split A: HF-anchored','Split B: LF-anchored','Location','best');

% Damping (same for both)
nexttile;
semilogx(Om, Ceq, 'LineWidth',1.8);
grid on; box on;
xline(Om_e,'k--');
xlabel('\Omega'); ylabel('C_{eq}(\Omega)');
title('Equivalent damping (unique)');

%% ---- Useful peak info ----
[Cmax, iC] = max(Ceq);
fprintf('Ceq peak = %.6f at Omega = %.6f\n', Cmax, Om(iC));