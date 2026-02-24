clc; clear; close all;
%求最优算子
%% ===================== User settings =====================
theta = 0.18;

% Frequency band
Om_min = 0.1;
Om_max = 10;
Nw = 350;
Om = logspace(log10(Om_min), log10(Om_max), Nw).';

% Split frequency (low/high band separator)
Om_split = 1.0;  % 你可以改成 0.8/1.5 等
idxL = (Om <= Om_split);
idxH = (Om >= Om_split);

% Search bounds for [sigma, kap_e, kap_c]
lb = [0.01, 0.01, 0.01];
ub = [2.00, 2.00, 2.00];

% Sampling settings
Nsamp = 4000;       % 先跑 2000~10000 都行（很快）
TopN  = 20;         % 输出前 TopN 组候选

% Trend targets (soft requirements)
% 你想：低频 Meq大, Ceq大, Keq小; 高频 Meq小, Ceq小, Keq大
% 这里用相对阈值：用候选自身的全频中位数作为参考，避免量纲问题
w = struct();
w.KL = 2.0;  % low-freq stiffness should be small
w.CL = 2.0;  % low-freq damping should be large
w.ML = 2.0;  % low-freq inertia should be large
w.KH = 2.0;  % high-freq stiffness should be large
w.CH = 1.5;  % high-freq damping should be small
w.MH = 1.5;  % high-freq inertia should be small

% Hard filters
minDenThresh = 0.03;      % avoid near-singularity
maxNegCproxyCount = 15;   % avoid too wide negative-damping proxy (optional)

rng(1); % reproducible

%% ===================== Sampling (LHS) =====================
X = lhsdesign(Nsamp,3);            % in (0,1)
X = lb + X.*(ub-lb);               % scale to bounds
% X(:,1)=sigma, X(:,2)=kap_e, X(:,3)=kap_c

J = inf(Nsamp,1);
feat = struct('CePeak',nan(Nsamp,1),'Om_e',nan(Nsamp,1), ...
              'KL',nan(Nsamp,1),'KH',nan(Nsamp,1), ...
              'CL',nan(Nsamp,1),'CH',nan(Nsamp,1), ...
              'ML',nan(Nsamp,1),'MH',nan(Nsamp,1));

%% ===================== Evaluate =====================
for i = 1:Nsamp
    sigma = X(i,1); kap_e = X(i,2); kap_c = X(i,3);

    % --- operator & denominator ---
    den = kap_e*Om.^2 - 1i*sigma*Om - kap_c;
    if min(abs(den)) < minDenThresh
        continue; % reject
    end

    % optional: negative damping proxy filter (very rough)
    K = (theta^2).*Om.^2 ./ den;
    Cproxy = imag(K)./max(Om,1e-12);   % equals Ceq actually (unique); will be >=0 if sigma>0
    if nnz(Cproxy < -0.2) > maxNegCproxyCount
        continue;
    end

    Kr  = real(K);
    Ceq = imag(K)./max(Om,1e-12);

    % --- LF-anchored split (recommended for your trend design) ---
    Meq = (theta^2)/kap_c * ones(size(Om));
    Keq = Kr + Om.^2 .* Meq;

    % --- compute band averages (log-weighted) ---
    % 用 log Ω 权重更符合你的对数频轴
    lw = ones(size(Om));
    lw = lw ./ sum(lw); % simple uniform weights on log grid (already log-spaced)

    meanL = @(y) sum(y(idxL).*lw(idxL)) / sum(lw(idxL));
    meanH = @(y) sum(y(idxH).*lw(idxH)) / sum(lw(idxH));

    KL = meanL(abs(Keq));
    KH = meanH(abs(Keq));
    CL = meanL(Ceq);
    CH = meanH(Ceq);
    ML = meanL(Meq);
    MH = meanH(Meq);

    % --- shape score design: smaller is better ---
    % 低频：K小 → penalize KL
    % 低频：C大 → penalize if CL is not large compared to CH
    % 低频：M大 → penalize if ML is not large compared to MH (here Meq constant so ratio=1)
    %
    % 高频：K大 → penalize if KH is not large compared to KL
    % 高频：C小 → penalize CH
    % 高频：M小 → penalize MH (constant, so it will push kap_c larger if you keep this term)
    %
    % 注意：Meq 在 LF-anchored 下是常数（频率不变），所以“低频大，高频小”无法同时满足。
    % 更合理的写法是：要求 Meq 大（用于低频惯性增强），同时不希望 KH 被它拖累——由 Keq 负责。
    %
    % 这里采取：鼓励 Meq 足够大（相对阈值），但不要求高频更小。
    Mtar = 1.0;  % 你可以改成你想要的惯容水平阈值
    Ctar = 0.05; % 阻尼低频平均的目标下限（按你自己的量级调）
    KtarH = 0.05; % 高频刚度平均的目标下限（量级按你系统调）

    % soft constraints via hinge loss
    hinge = @(z) max(0,z);

    Jshape = 0;
    Jshape = Jshape + w.KL * KL;
    Jshape = Jshape + w.CL * hinge(Ctar - CL);
    Jshape = Jshape + w.ML * hinge(Mtar - ML);

    Jshape = Jshape + w.KH * hinge(KtarH - KH);
    Jshape = Jshape + w.CH * CH;
    % 不再强迫 MH 小（否则会和你“低频惯容大”矛盾）
    % Jshape = Jshape + w.MH * MH;

    % add a mild term to prefer wider damping enhancement:
    % penalize overly sharp Ceq peak: ratio peak/meanL too large
    [CePeak,~] = max(Ceq);
    sharp = CePeak / max(CL,1e-12);
    Jshape = Jshape + 0.3 * hinge(sharp - 15);

    J(i) = Jshape;

    % store features
    feat.CePeak(i) = CePeak;
    feat.Om_e(i)   = sqrt(kap_c/kap_e);
    feat.KL(i) = KL; feat.KH(i) = KH;
    feat.CL(i) = CL; feat.CH(i) = CH;
    feat.ML(i) = ML; feat.MH(i) = MH;
end

%% ===================== Rank & report =====================
[Js, idx] = sort(J, 'ascend');
idx = idx(isfinite(Js));
idx = idx(1:min(TopN,numel(idx)));

fprintf('\n==== Top-%d candidates (operator-domain) ====\n', numel(idx));
fprintf(' rank |     Jshape |  sigma   kap_e   kap_c |  Om_e  |  Ce_peak  |   KL     KH     CL     CH     Meq\n');
fprintf('--------------------------------------------------------------------------------------------------------\n');

for r = 1:numel(idx)
    i = idx(r);
    sigma = X(i,1); kap_e = X(i,2); kap_c = X(i,3);
    fprintf('%4d | %9.4e | %6.3f  %6.3f  %6.3f | %5.3f | %8.4f | %6.3f %6.3f %6.3f %6.3f %6.3f\n', ...
        r, J(i), sigma, kap_e, kap_c, feat.Om_e(i), feat.CePeak(i), ...
        feat.KL(i), feat.KH(i), feat.CL(i), feat.CH(i), feat.ML(i));
end

%% ===================== Plot best candidate curves =====================
best = idx(1);
sigma = X(best,1); kap_e = X(best,2); kap_c = X(best,3);

den = kap_e*Om.^2 - 1i*sigma*Om - kap_c;
K   = (theta^2).*Om.^2 ./ den;
Kr  = real(K);
Ki  = imag(K);
Ceq = Ki ./ max(Om,1e-12);

Meq = (theta^2)/kap_c * ones(size(Om));
Keq = Kr + Om.^2 .* Meq;

Om_e = sqrt(kap_c/kap_e);

figure('Color','w','Position',[120 120 1100 420]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile;
semilogx(Om, Meq,'LineWidth',1.6); grid on; box on;
xline(Om_e,'k--','\Omega_e'); xline(Om_split,'k:','\Omega_s');
xlabel('\Omega'); ylabel('M_{eq}(\Omega)');
title(sprintf('Meq (LF-anchored), \\sigma=%.3f, \\kappa_e=%.3f, \\kappa_c=%.3f',sigma,kap_e,kap_c));

nexttile;
semilogx(Om, Keq,'LineWidth',1.6); grid on; box on;
xline(Om_e,'k--','\Omega_e'); xline(Om_split,'k:','\Omega_s');
xlabel('\Omega'); ylabel('K_{eq}(\Omega)');
title('Keq');

nexttile;
semilogx(Om, Ceq,'LineWidth',1.6); grid on; box on;
xline(Om_e,'k--','\Omega_e'); xline(Om_split,'k:','\Omega_s');
xlabel('\Omega'); ylabel('C_{eq}(\Omega)');
title('Ceq (unique)');

figure('Color','w','Position',[150 150 900 380]);
semilogx(Om, Kr,'LineWidth',1.5); hold on;
semilogx(Om, Ki,'LineWidth',1.5);
grid on; box on;
xline(Om_e,'k--','\Omega_e'); xline(Om_split,'k:','\Omega_s');
xlabel('\Omega'); ylabel('\mathcal{K}_r,\ \mathcal{K}_i');
title('\mathcal{K}(\Omega) real/imag');
legend('\mathcal{K}_r','\mathcal{K}_i','Location','best');