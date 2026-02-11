%% Run_Sweep_kapc.m
% 参数敏感性分析 - RLC 中 kap_c (sysP(9))
% 目的：
% 1) 对比 kap_c 对频响曲线的影响
% 2) 纵坐标使用 dB (20log10)
% 3) 自动输出两个主要共振峰（峰值/峰位）
%
% 依赖：FRF.m, nondim_temp2.m, newton.m, branch_follow2.m 等在同一路径

clear; clc; close all;

global Fw FixedOmega
Fw = 0.05;
FixedOmega = [];

%% 1) 基础参数（除 kap_c 外固定）
P.be1   = 1.0;   % sysP(1)
P.be2   = 1.0;   % sysP(2)
P.mu    = 0.2;   % sysP(3)
P.al1   = 0.01;   % sysP(4)
P.ga1   = 0.1;   % sysP(5)
P.ze1   = 0.05;  % sysP(6)
P.lam   = 0.5

% 调谐后的电路参数（按你当前方案）
P.kap_e = 2.04;  % sysP(8)
%P.kap_c = 1.0;   % sysP(9)
P.sigma = 0.5;   % sysP(10)
P.ga2   = 0.2;   % sysP(11)
%% 2) kap_c 扫描列表（建议跨数量级）
kapc_list = [0.2, 0.5, 1.0, 2.0, 5];

styles = {'k-', 'b-', 'g-', 'r-', 'm-'};
mkStep = 20;

%% 3) 画图
figure('Position',[100,100,1100,450],'Color','w');

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
xlabel(ax1, '\Omega'); ylabel(ax1, '20log_{10}(|X_1|) (dB)');
title(ax1, 'Upper Mass FRF vs \kappa_c');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
xlabel(ax2, '\Omega'); ylabel(ax2, '20log_{10}(|Q|) (dB)');
title(ax2, 'Circuit State |Q| vs \kappa_c');

legend_str = cell(numel(kapc_list),1);

fprintf('================================================\n');
fprintf('SWEEP: kap_c = %s\n', mat2str(kapc_list));
fprintf('Fw = %.4f, lambda = %.4f\n', Fw, P.lam);
fprintf('kap_e = %.4f, sigma = %.4f\n', P.kap_e, P.sigma);
fprintf('================================================\n');

PeakTable = [];  % [kap_c, mode_id, Omega_peak, AmpX1_dB_peak]

for k = 1:numel(kapc_list)
    P.kap_c = kapc_list(k);

    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
            P.kap_e, P.kap_c, P.sigma, P.ga2];

    fprintf('\n>>> Case %d/%d: kap_c = %.4f\n', k, numel(kapc_list), P.kap_c);

    try
        x_res = FRF(sysP);   % 16xN, last row Omega
    catch ME
        fprintf('   [ERROR] FRF failed: %s\n', ME.message);
        continue;
    end

    if isempty(x_res) || size(x_res,1) < 16
        fprintf('   [WARN] Empty/invalid FRF result.\n');
        continue;
    end

    Omega = x_res(16, :);

    AmpX1 = sqrt(x_res(2,:).^2  + x_res(3,:).^2);
    AmpQ  = sqrt(x_res(12,:).^2 + x_res(13,:).^2);

    AmpX1_dB = 20*log10(AmpX1 + eps);
    AmpQ_dB  = 20*log10(AmpQ  + eps);

    p1 = plot(ax1, Omega, AmpX1_dB, styles{k}, 'LineWidth', 1.8);
    p2 = plot(ax2, Omega, AmpQ_dB,  styles{k}, 'LineWidth', 1.8);

    if numel(Omega) > mkStep
        p1.Marker = '.';
        p1.MarkerIndices = 1:mkStep:numel(Omega);
        p2.Marker = '.';
        p2.MarkerIndices = 1:mkStep:numel(Omega);
    end

    legend_str{k} = sprintf('\\kappa_c = %.3g', P.kap_c);

    % 两个主峰（更稳：用 findpeaks）
    [pk_dB, pk_w] = find_two_main_peaks_findpeaks(AmpX1_dB, Omega);

    if ~isempty(pk_dB)
        for m = 1:numel(pk_dB)
            PeakTable = [PeakTable; P.kap_c, m, pk_w(m), pk_dB(m)]; %#ok<AGROW>
            fprintf('   Peak %d: Omega = %.4f, |X1| = %.2f dB\n', m, pk_w(m), pk_dB(m));
        end
    else
        fprintf('   [WARN] Peaks not found.\n');
    end
end

legend(ax1, legend_str, 'Location','best');
legend(ax2, legend_str, 'Location','best');

% 不强行截断，避免漏掉峰
xlim(ax1, [min(Omega) max(Omega)]);
xlim(ax2, [min(Omega) max(Omega)]);

fprintf('\n================= Peak Summary (|X1| dB) =================\n');
fprintf(' kap_c   mode     Omega_peak      Amp_dB\n');
fprintf('-----------------------------------------------------------\n');
for i = 1:size(PeakTable,1)
    fprintf(' %6.3f   %4d     %10.4f    %8.2f\n', PeakTable(i,1), PeakTable(i,2), PeakTable(i,3), PeakTable(i,4));
end
fprintf('===========================================================\n');

%% -------- 辅助：findpeaks版更稳 --------
function [pk_dB, pk_w] = find_two_main_peaks_findpeaks(y_dB, w)
    pk_dB = []; pk_w = [];
    if numel(y_dB) < 20, return; end

    try
        [pks, locs] = findpeaks(y_dB, w, ...
            'MinPeakProminence', 2, ...   % 2~5 dB 可调
            'MinPeakDistance',   0.2);
    catch
        return;
    end
    if isempty(pks), return; end

    [~, ord] = sort(pks, 'descend');
    ord = ord(1:min(2,numel(ord)));

    [locs2, ord2] = sort(locs(ord), 'ascend');
    pk_w  = locs2(:).';
    pk_dB = pks(ord(ord2)).';
end
