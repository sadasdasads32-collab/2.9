%% Run_Sweep_Lambda_Dissipation.m
% 目的：
%  1) 扫描 lambda，计算并绘制电阻耗散功率指标 <P_R>=sigma*<Q'^2>
%  2) 与机械通道(X1/X2/Xrel)的 FRF 同步对照
%  3) 自动输出每个 lambda 的最大耗散值及其对应频率
%
% 依赖：FRF.m, nondim_temp2.m, newton.m, branch_follow2.m 等

clear; clc; close all;

global Fw FixedOmega
Fw = 0.05;
FixedOmega = [];

%% 1) 固定参数（除 lam）
P.be1   = 1.0;
P.be2   = 1.0;
P.mu    = 0.2;      % 你可改回 0.2 / 0.8 对比
P.al1   = 0.0;
P.ga1   = 0.1;
P.ze1   = 0.05;
P.kap_e = 0.05;
P.kap_c = 1.0;
P.sigma = 2.0;      % sigma>0 才是“电阻耗散”
P.ga2   = 0.2;

%% 2) lambda 列表
lam_list = [0, 0.2, 0.3, 0.4, 0.5];
styles   = {'k-', 'b-', 'g-', 'r-', 'm-'};

%% 3) 频段设置
Range = [0.2, 4.5];          % 全频段
RangeHF = [1.2, 2.4];        % 你也可以指定“主作用频段”用于打印/放大

fprintf('================================================\n');
fprintf('Fw=%.4f, mu=%.3f, sigma=%.3f, lambda sweep=%s\n', Fw, P.mu, P.sigma, mat2str(lam_list));
fprintf('================================================\n');

%% 4) 图：机械(三通道) + 电路耗散(一张) + 可选|Q|
figure('Position',[80,60,1000,900],'Color','w');
ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
ylabel(ax1,'20log_{10}(|X_1|) (dB)'); title(ax1,'|X_1| vs \lambda');

ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
ylabel(ax2,'20log_{10}(|X_2|) (dB)'); title(ax2,'|X_2| vs \lambda');

ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
xlabel(ax3,'\Omega'); ylabel(ax3,'20log_{10}(|X_{rel}|) (dB)'); title(ax3,'|X_{rel}=X_1-X_2| vs \lambda');

figure('Position',[1120,120,720,520],'Color','w');
axD = axes; hold(axD,'on'); grid(axD,'on'); box(axD,'on');
xlabel(axD,'\Omega'); ylabel(axD,'10log_{10}(\langle P_R \rangle)  (dB)');  % 功率用10log
title(axD,'Resistor dissipation:  \langle P_R \rangle = \sigma \cdot \Omega^2 (Q_c^2+Q_s^2)/2');

figure('Position',[1120,680,720,420],'Color','w');
axQ = axes; hold(axQ,'on'); grid(axQ,'on'); box(axQ,'on');
xlabel(axQ,'\Omega'); ylabel(axQ,'20log_{10}(|Q|) (dB)');
title(axQ,'Circuit state |Q| vs \lambda');

leg = strings(numel(lam_list),1);

%% 5) 汇总表：lambda, Omega_at_maxD, maxD, maxD_dB, maxD_in_RangeHF
Summary = [];

for k = 1:numel(lam_list)
    P.lam = lam_list(k);

    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
            P.kap_e, P.kap_c, P.sigma, P.ga2];

    fprintf('\n>>> lambda = %.2f\n', P.lam);

    % --- FRF ---
    try
        x_res = FRF(sysP);
    catch ME
        fprintf('   [ERROR] FRF failed: %s\n', ME.message);
        continue;
    end

    if isempty(x_res) || size(x_res,1) < 16 || size(x_res,2) < 10
        fprintf('   [WARN] FRF output invalid.\n');
        continue;
    end

    Omega = x_res(16,:);

    % --- 基波系数（与你之前一致） ---
    x1c = x_res(2,:);  x1s = x_res(3,:);
    x2c = x_res(7,:);  x2s = x_res(8,:);
    Qc1 = x_res(12,:); Qs1 = x_res(13,:);

    % --- 幅值(dB) ---
    X1_dB   = 20*log10(hypot(x1c, x1s) + eps);
    X2_dB   = 20*log10(hypot(x2c, x2s) + eps);
    Xrel_dB = 20*log10(hypot(x1c-x2c, x1s-x2s) + eps);
    Q_dB    = 20*log10(hypot(Qc1, Qs1) + eps);

    % --- 电阻耗散：平均功率 <P_R> = sigma * <Q'^2> ---
    % Q' 的基波 RMS^2 = Omega^2*(Qc^2+Qs^2)/2
    PR_avg = P.sigma * (Omega.^2) .* (Qc1.^2 + Qs1.^2) / 2;

    % 功率dB：用10log10
    PR_dB = 10*log10(PR_avg + eps);

    % --- 画图 ---
    plot(ax1, Omega, X1_dB, styles{k}, 'LineWidth', 1.8);
    plot(ax2, Omega, X2_dB, styles{k}, 'LineWidth', 1.8);
    plot(ax3, Omega, Xrel_dB, styles{k}, 'LineWidth', 1.8);

    plot(axD, Omega, PR_dB, styles{k}, 'LineWidth', 1.8);
    plot(axQ, Omega, Q_dB,  styles{k}, 'LineWidth', 1.8);

    leg(k) = sprintf('\\lambda=%.1f', P.lam);

    % --- 输出“最大耗散”位置（全频）---
    [maxD, idx] = max(PR_avg);
    Om_maxD = Omega(idx);
    maxD_dB = 10*log10(maxD + eps);

    % --- 输出“指定频段RangeHF内最大耗散”（更有针对性）---
    maskHF = (Omega>=RangeHF(1)) & (Omega<=RangeHF(2));
    if any(maskHF)
        [maxD_HF, idxHF_local] = max(PR_avg(maskHF));
        OmHF = Omega(maskHF);
        Om_maxD_HF = OmHF(idxHF_local);
        maxD_HF_dB = 10*log10(maxD_HF + eps);
    else
        Om_maxD_HF = NaN; maxD_HF = NaN; maxD_HF_dB = NaN;
    end

    Summary = [Summary; P.lam, Om_maxD, maxD, maxD_dB, Om_maxD_HF, maxD_HF, maxD_HF_dB]; %#ok<AGROW>

    fprintf('   Max dissipation (all):   Omega=%.4f, <PR>=%.3e (%.2f dB)\n', Om_maxD, maxD, maxD_dB);
    fprintf('   Max dissipation (%.1f-%.1f): Omega=%.4f, <PR>=%.3e (%.2f dB)\n', ...
        RangeHF(1), RangeHF(2), Om_maxD_HF, maxD_HF, maxD_HF_dB);
end

% 统一横轴
xlim(ax1, Range); xlim(ax2, Range); xlim(ax3, Range);
xlim(axD, Range); xlim(axQ, Range);

legend(ax1, leg, 'Location','best');
legend(ax2, leg, 'Location','best');
legend(ax3, leg, 'Location','best');
legend(axD, leg, 'Location','best');
legend(axQ, leg, 'Location','best');

%% 6) 汇总表打印
fprintf('\n================ Dissipation Summary =================\n');
fprintf('lambda   Om_max(all)   <PR>_max(all)     dB    Om_max(HF)   <PR>_max(HF)      dB\n');
fprintf('---------------------------------------------------------------------------------\n');
for i=1:size(Summary,1)
    fprintf('%5.2f   %10.4f   %12.3e   %7.2f   %10.4f   %12.3e   %7.2f\n', ...
        Summary(i,1), Summary(i,2), Summary(i,3), Summary(i,4), Summary(i,5), Summary(i,6), Summary(i,7));
end
fprintf('======================================================\n');

% 额外提醒：sigma>0 时 PR_avg 应>=0（数值误差除外）
if P.sigma > 0
    fprintf('\n[Check] sigma>0: <PR> should be non-negative. If you see negative values, check indices or sigma sign.\n');
end
