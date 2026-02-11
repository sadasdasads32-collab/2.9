%% Run_Sweep_Lambda.m
% 参数敏感性分析 - 机电耦合强度 Lambda (lam)
% 目的：
% 1) 对比 lam = [0, ...] 对频响曲线的影响
% 2) 纵坐标使用 dB (20log10)
% 3) 自动输出两个主要共振峰（峰值/峰位）
% 4) 同时输出电路“能量/耗散通道”指标：|Q'| ≈ |\Omega Q|
%    并给出耗散代理指标 sigma*max(|Q'|^2)
%
% 输出：
%  - Figure: 左 |X1|dB vs Omega，右 |\Omega Q|dB vs Omega（实线）
%           同色虚线为 |Q|dB 作为参考
%  - PeakTable: [lambda, mode_id, Omega_peak, AmpX1_dB_peak]
%  - DissTable: [lambda, sigma*max(|Q'|^2)]
%
% 依赖：FRF.m, nondim_temp2.m, newton.m, branch_follow2.m 等在同一路径

clear; clc; close all;

%% -----------------------------
% 0) 全局变量（与你的 nondim_temp2 约定一致）
% -----------------------------
global Fw FixedOmega
Fw = 0.05;        % 激励幅值（保持不变，方便对比）
FixedOmega = [];  % 扫频模式

%% -----------------------------
% 1) 基础参数（除 lam 外固定）
% -----------------------------
P.be1   = 1.0;   % sysP(1)
P.be2   = 1.0;   % sysP(2)
P.mu    = 0.2;   % sysP(3)
P.al1   = 0.01;   % sysP(4)
P.ga1   = 0.1;   % sysP(5)
P.ze1   = 0.05;  % sysP(6)
% lam 是待扫描变量 sysP(7)

% 调谐后的电路参数（按你当前方案）
P.kap_e = 2.04;  % sysP(8)
P.kap_c = 1.0;   % sysP(9)
P.sigma = 0.5;   % sysP(10)
P.ga2   = 0.2;   % sysP(11)

%% -----------------------------
% 2) lam 扫描列表
% -----------------------------
lam_list = [0, 0.2, 0.3, 0.4, 0.5];

% 线型/颜色
styles = {'k-', 'b-', 'g-', 'r-', 'm-'};
mkStep = 20;

%% -----------------------------
% 3) 画图
% -----------------------------
figure('Position',[100,100,1100,450],'Color','w');

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
xlabel(ax1, '\Omega'); ylabel(ax1, '20log_{10}(|X_1|)  (dB)');
title(ax1, 'Upper Mass FRF vs \lambda');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
xlabel(ax2, '\Omega'); ylabel(ax2, '20log_{10}(|\Omega Q|)  (dB)');
title(ax2, 'Circuit Energy Proxy |\Omega Q| (solid) & |Q| (dashed) vs \lambda');

legend_str = cell(numel(lam_list),1);

fprintf('================================================\n');
fprintf('SWEEP: lambda = %s\n', mat2str(lam_list));
fprintf('Fw = %.4f\n', Fw);
fprintf('Fixed params: kap_e=%.4f, kap_c=%.4f, sigma=%.4f\n', P.kap_e, P.kap_c, P.sigma);
fprintf('================================================\n');

% 汇总表
PeakTable = [];  % [lam, mode_id, Omega_peak, AmpX1_dB_peak]
DissTable = [];  % [lam, sigma*max(|Q'|^2)]  (Q'≈Omega*Q)

Omega_last = [];

%% -----------------------------
% 4) 扫描循环
% -----------------------------
for k = 1:numel(lam_list)
    P.lam = lam_list(k);

    % sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
            P.kap_e, P.kap_c, P.sigma, P.ga2];

    fprintf('\n>>> Case %d/%d: lambda = %.2f\n', k, numel(lam_list), P.lam);

    % --- FRF 扫频 ---
    try
        x_res = FRF(sysP);   % 16xN，最后一行为 Omega
    catch ME
        fprintf('   [ERROR] FRF failed: %s\n', ME.message);
        continue;
    end

    if isempty(x_res) || size(x_res,1) < 16
        fprintf('   [WARN] Empty/invalid FRF result.\n');
        continue;
    end

    Omega = x_res(16, :);
    Omega_last = Omega;

    % 幅值（1次谐波）：x = a*cos + b*sin => Amp = sqrt(a^2 + b^2)
    AmpX1 = sqrt(x_res(2,:).^2  + x_res(3,:).^2);
    AmpQ  = sqrt(x_res(12,:).^2 + x_res(13,:).^2);

    % dB
    AmpX1_dB = 20*log10(AmpX1 + eps);
    AmpQ_dB  = 20*log10(AmpQ  + eps);

    % 能量/耗散通道：|Q'| ≈ |\Omega Q|
    AmpQp    = abs(Omega) .* AmpQ;           % |Q'| 幅值近似
    AmpQp_dB = 20*log10(AmpQp + eps);

    % 耗散代理指标：sigma*max(|Q'|^2)
    DissMetric = P.sigma * max(AmpQp.^2);
    DissTable = [DissTable; P.lam, DissMetric]; %#ok<AGROW>
    fprintf('   Dissipation proxy: sigma*max(|Q''|^2) = %.3e\n', DissMetric);

    % --- 画图：左 X1，右 OmegaQ(实线) + Q(虚线) ---
    p1 = plot(ax1, Omega, AmpX1_dB, styles{k}, 'LineWidth', 1.8);
    p2 = plot(ax2, Omega, AmpQp_dB, styles{k}, 'LineWidth', 1.8); % 实线：|\Omega Q|
    plot(ax2, Omega, AmpQ_dB,  [styles{k}(1) '--'], 'LineWidth', 1.2); % 虚线：|Q|

    % 打点
    if numel(Omega) > mkStep
        p1.Marker = '.';
        p1.MarkerIndices = 1:mkStep:numel(Omega);

        p2.Marker = '.';
        p2.MarkerIndices = 1:mkStep:numel(Omega);
    end

    legend_str{k} = sprintf('\\lambda = %.1f', P.lam);

    % --- 寻峰：找两个主峰（findpeaks 稳健版） ---
    [pk_dB, pk_w] = find_two_main_peaks_findpeaks(AmpX1_dB, Omega);

    if ~isempty(pk_dB)
        for m = 1:numel(pk_dB)
            PeakTable = [PeakTable; P.lam, m, pk_w(m), pk_dB(m)]; %#ok<AGROW>
            fprintf('   Peak %d: Omega = %.4f, |X1| = %.2f dB\n', m, pk_w(m), pk_dB(m));
        end
    else
        fprintf('   [WARN] Peaks not found.\n');
    end
end

%% -----------------------------
% 5) 图例、坐标范围
% -----------------------------
legend(ax1, legend_str, 'Location','best');
legend(ax2, legend_str, 'Location','best');

% 你可以固定范围，也可以自动范围（推荐自动避免漏峰）
if ~isempty(Omega_last)
    xlim(ax1, [0.2, 4.5]);
    xlim(ax2, [0.2, 4.5]);
end

%% -----------------------------
% 6) 打印表格
% -----------------------------
fprintf('\n================= Peak Summary (|X1| dB) =================\n');
fprintf(' lambda   mode     Omega_peak      Amp_dB\n');
fprintf('-----------------------------------------------------------\n');
for i = 1:size(PeakTable,1)
    fprintf(' %6.2f   %4d     %10.4f    %8.2f\n', PeakTable(i,1), PeakTable(i,2), PeakTable(i,3), PeakTable(i,4));
end
fprintf('===========================================================\n');

fprintf('\n============= Dissipation Proxy Summary =============\n');
fprintf(' lambda     sigma*max(|Q''|^2)\n');
fprintf('-----------------------------------------------------\n');
for i = 1:size(DissTable,1)
    fprintf(' %6.2f     %12.4e\n', DissTable(i,1), DissTable(i,2));
end
fprintf('=====================================================\n');

%% -----------------------------
% 辅助函数：找两个主峰（findpeaks 稳健版）
% -----------------------------
function [pk_dB, pk_w] = find_two_main_peaks_findpeaks(y_dB, w)
    pk_dB = []; pk_w = [];
    if numel(y_dB) < 20, return; end

    % 参数可按你的曲线陡峭程度调整
    MinProm = 2;     % 最小“显著性”(dB)
    MinDist = 0.15;  % 最小峰距(频率)

    try
        [pks, locs] = findpeaks(y_dB, w, ...
            'MinPeakProminence', MinProm, ...
            'MinPeakDistance',   MinDist);
    catch
        return;
    end

    if isempty(pks), return; end

    % 取最大两个峰
    [~, ord] = sort(pks, 'descend');
    ord = ord(1:min(2, numel(ord)));

    % 再按频率升序（Mode1/Mode2）
    [pk_w, idx2] = sort(locs(ord), 'ascend');
    pk_dB = pks(ord(idx2));
end
