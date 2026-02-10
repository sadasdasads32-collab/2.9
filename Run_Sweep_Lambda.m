%% Run_Sweep_Lambda.m
% 参数敏感性分析 - 机电耦合强度 Lambda (lam)
% 目的：
% 1) 对比 lam = [0, 0.2, 0.4, 0.6, 0.8] 对频响曲线的影响
% 2) 纵坐标使用 dB (20log10)
% 3) 自动输出两个主要共振峰（峰值/峰位）
%
% 依赖：FRF.m, nondim_temp2.m, newton.m, branch_follow2.m 等在同一路径

clear; clc; close all;

% -----------------------------
% 0) 全局变量（与你的 nondim_temp2 约定一致）
% -----------------------------
global Fw FixedOmega
Fw = 0.05;        % 激励幅值（保持不变，方便对比）
FixedOmega = [];  % 扫频模式

% -----------------------------
% 1) 基础参数（除 lam 外全部固定）
%    你可以按你论文参数替换这里
% -----------------------------
P.be1   = 1.0;   % sysP(1)
P.be2   = 1.0;   % sysP(2)
P.mu    = 0.2;   % sysP(3) - 修正为文献常用值
P.al1   = 0.0;   % sysP(4) - 上层线性刚度
P.ga1   = 0.1;   % sysP(5) - 上层非线性刚度
P.ze1   = 0.05;  % sysP(6) - 下层阻尼比
% lam 是待扫描变量 sysP(7)
P.kap_e = 0.05;  % sysP(8)
P.kap_c = 1.0;   % sysP(9)
P.sigma = 2.0;   % sysP(10)- 分流电阻 
P.ga2   = 0.2;   % sysP(11)- 下层非线性


% -----------------------------
% 2) lam 扫描列表
% -----------------------------
lam_list = [0, 0.2, 0.3, 0.4, 0.5];

% 颜色/线型（保证可区分）
styles = {'k-', 'b-', 'g-', 'r-', 'm-'};
mkStep = 20; % 每隔多少点打一个 marker（可调）

% -----------------------------
% 3) 画图：|X1| dB vs Omega
% -----------------------------
figure('Position',[100,100,1100,450],'Color','w');

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
xlabel(ax1, '\Omega'); ylabel(ax1, '20log_{10}(|X_1|)  (dB)');
title(ax1, 'Upper Mass FRF vs \lambda');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
xlabel(ax2, '\Omega'); ylabel(ax2, '20log_{10}(|Q|)  (dB)');
title(ax2, 'Circuit State |Q| vs \lambda');

legend_str = cell(numel(lam_list),1);

fprintf('================================================\n');
fprintf('SWEEP: lambda = %s\n', mat2str(lam_list));
fprintf('Fw = %.4f\n', Fw);
fprintf('================================================\n');

% 用于汇总打印（每个 lambda 两个主峰）
PeakTable = [];  % columns: [lam, mode_id, Omega_peak, AmpX1_dB_peak]

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

    % 幅值（1次谐波）：x = a*cos + b*sin => Amp = sqrt(a^2 + b^2)
    AmpX1 = sqrt(x_res(2,:).^2  + x_res(3,:).^2);
    AmpQ  = sqrt(x_res(12,:).^2 + x_res(13,:).^2);

    % 转 dB（加 eps 防止 log(0)）
    AmpX1_dB = 20*log10(AmpX1 + eps);
    AmpQ_dB  = 20*log10(AmpQ  + eps);

    % --- 画图 ---
    p1 = plot(ax1, Omega, AmpX1_dB, styles{k}, 'LineWidth', 1.8);
    p2 = plot(ax2, Omega, AmpQ_dB,  styles{k}, 'LineWidth', 1.8);

    % 可选：打点
    if numel(Omega) > mkStep
        p1.Marker = '.';
        p1.MarkerIndices = 1:mkStep:numel(Omega);

        p2.Marker = '.';
        p2.MarkerIndices = 1:mkStep:numel(Omega);
    end

    legend_str{k} = sprintf('\\lambda = %.1f', P.lam);

    % --- 寻峰：只找“最主要的两个峰”（通常对应 1、2 阶模态） ---
    [pk_dB, pk_w] = find_two_main_peaks(AmpX1_dB, Omega);

    if ~isempty(pk_dB)
        for m = 1:numel(pk_dB)
            PeakTable = [PeakTable; P.lam, m, pk_w(m), pk_dB(m)];
            fprintf('   Peak %d: Omega = %.4f, |X1| = %.2f dB\n', m, pk_w(m), pk_dB(m));
        end
    else
        fprintf('   [WARN] Peaks not found (maybe curve monotonic or too few points).\n');
    end
end

legend(ax1, legend_str, 'Location','best');
legend(ax2, legend_str, 'Location','best');

% 统一横轴范围（按你扫频结果可调）
xlim(ax1, [0.2, 4.5]);
xlim(ax2, [0.2, 4.5]);

% 打印汇总表
fprintf('\n================= Peak Summary (|X1| dB) =================\n');
fprintf(' lambda   mode     Omega_peak      Amp_dB\n');
fprintf('-----------------------------------------------------------\n');
for i = 1:size(PeakTable,1)
    fprintf(' %6.2f   %4d     %10.4f    %8.2f\n', PeakTable(i,1), PeakTable(i,2), PeakTable(i,3), PeakTable(i,4));
end
fprintf('===========================================================\n');

% -----------------------------
% 辅助函数：找两个主峰（简单稳健版）
% -----------------------------
function [pk_dB, pk_w] = find_two_main_peaks(y_dB, w)
    pk_dB = []; pk_w = [];
    if numel(y_dB) < 5, return; end

    % 只找局部极大值
    idx = [];
    for i = 2:numel(y_dB)-1
        if y_dB(i) > y_dB(i-1) && y_dB(i) > y_dB(i+1)
            idx(end+1) = i; %#ok<AGROW>
        end
    end
    if isempty(idx), return; end

    % 把峰按峰值从大到小排序，取前 2 个
    [~, ord] = sort(y_dB(idx), 'descend');
    idx = idx(ord);

    take = min(2, numel(idx));
    idx = idx(1:take);

    % 再按频率从小到大排列（对应 Mode1/Mode2）
    [~, ord2] = sort(w(idx), 'ascend');
    idx = idx(ord2);

    pk_dB = y_dB(idx);
    pk_w  = w(idx);
end
