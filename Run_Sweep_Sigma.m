%% Run_Sweep_Sigma.m
% 参数敏感性分析 - RLC 中分流电阻 Sigma (sigma)
% 目的：
% 1) 对比 sigma = [..] 对频响曲线的影响
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
% 1) 基础参数（除 sigma 外全部固定）
P.be1   = 1.0;   % sysP(1)
P.be2   = 1.0;   % sysP(2)
P.mu    = 0.2;   % sysP(3)
P.al1   = 0.01;   % sysP(4)
P.ga1   = 0.1;   % sysP(5)
P.ze1   = 0.05;  % sysP(6)
P.lam   = 0.5

% 调谐后的电路参数（按你当前方案）
P.kap_e = 2.04;  % sysP(8)
P.kap_c = 1.0;   % sysP(9)
%P.sigma = 0.5;   % sysP(10)
P.ga2   = 0.2;   % sysP(11)

% -----------------------------
% 2) sigma 扫描列表（你可以改得更密）
% -----------------------------
sigma_list = [0.1, 0.5, 1.0, 2.0, 5.0];  % 推荐：跨数量级更容易看趋势
styles = {'k-', 'b-', 'g-', 'r-', 'm-'};
mkStep = 20; % 每隔多少点打一个 marker（可调）

% -----------------------------
% 3) 画图：|X1| dB vs Omega
% -----------------------------
figure('Position',[100,100,1100,450],'Color','w');

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
xlabel(ax1, '\Omega'); ylabel(ax1, '20log_{10}(|X_1|)  (dB)');
title(ax1, 'Upper Mass FRF vs \sigma');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
xlabel(ax2, '\Omega'); ylabel(ax2, '20log_{10}(|Q|)  (dB)');
title(ax2, 'Circuit State |Q| vs \sigma');

legend_str = cell(numel(sigma_list),1);

fprintf('================================================\n');
fprintf('SWEEP: sigma = %s\n', mat2str(sigma_list));
fprintf('Fw = %.4f, lambda = %.4f\n', Fw, P.lam);
fprintf('kap_e = %.4f, kap_c = %.4f\n', P.kap_e, P.kap_c);
fprintf('================================================\n');

% 汇总表：每个 sigma 两个主峰
PeakTable = [];  % columns: [sigma, mode_id, Omega_peak, AmpX1_dB_peak]

for k = 1:numel(sigma_list)
    P.sigma = sigma_list(k);

    % sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
    sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, ...
            P.kap_e, P.kap_c, P.sigma, P.ga2];

    fprintf('\n>>> Case %d/%d: sigma = %.4f\n', k, numel(sigma_list), P.sigma);

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

    if numel(Omega) > mkStep
        p1.Marker = '.';
        p1.MarkerIndices = 1:mkStep:numel(Omega);
        p2.Marker = '.';
        p2.MarkerIndices = 1:mkStep:numel(Omega);
    end

    legend_str{k} = sprintf('\\sigma = %.3g', P.sigma);

    % --- 寻峰：找两个主峰 ---
    [pk_dB, pk_w] = find_two_main_peaks(AmpX1_dB, Omega);

    if ~isempty(pk_dB)
        for m = 1:numel(pk_dB)
            PeakTable = [PeakTable; P.sigma, m, pk_w(m), pk_dB(m)]; %#ok<AGROW>
            fprintf('   Peak %d: Omega = %.4f, |X1| = %.2f dB\n', m, pk_w(m), pk_dB(m));
        end
    else
        fprintf('   [WARN] Peaks not found (maybe curve monotonic or too few points).\n');
    end
    % ---- 保存参考曲线（比如 sigma=2 时）----
    if abs(P.sigma - 2.0) < 1e-12
        AmpX1_dB_ref = AmpX1_dB;
        Omega_ref = Omega;
    end
    
    % ---- 在拿到参考后，比较差异 ----
    if exist('AmpX1_dB_ref','var')
        % 如果 Omega 网格一致，直接减；不一致就先插值
        if numel(Omega)==numel(Omega_ref) && max(abs(Omega-Omega_ref))<1e-10
            Delta = AmpX1_dB - AmpX1_dB_ref;
        else
            Delta = AmpX1_dB - interp1(Omega_ref, AmpX1_dB_ref, Omega, 'linear','extrap');
        end
        fprintf('   max |ΔX1| vs sigma=2 : %.2f dB\n', max(abs(Delta)));
    end

end

legend(ax1, legend_str, 'Location','best');
legend(ax2, legend_str, 'Location','best');

% 统一横轴范围（按你扫频结果可调）
xlim(ax1, [min(Omega) max(Omega)]);
xlim(ax2, [min(Omega) max(Omega)]);


% 打印汇总表
fprintf('\n================= Peak Summary (|X1| dB) =================\n');
fprintf(' sigma    mode     Omega_peak      Amp_dB\n');
fprintf('-----------------------------------------------------------\n');
for i = 1:size(PeakTable,1)
    fprintf(' %7.4f   %4d     %10.4f    %8.2f\n', PeakTable(i,1), PeakTable(i,2), PeakTable(i,3), PeakTable(i,4));
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
