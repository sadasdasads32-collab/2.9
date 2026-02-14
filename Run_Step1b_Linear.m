%% Run_Step1b_Linear_Verification.m
% 审稿人特供：第二步验证 - 线性纯机械系统 (dB版 + 自动寻峰 + 理论固有频率对照)
% 目的：
% 1) 验证线性模态 (Linear Modes) 是否合理
% 2) 打印共振峰坐标，量化固有频率
% 3) 与经典 two-stage（2DOF）线性模型理论固有频率对照

clear; clc; close all;

%% -----------------------------
% 1) 参数设置：线性纯机械（关键：电路参数全部置零）
% -----------------------------
P.be1 = 1;  %上层刚度与（基准刚度也就是上层刚度之比）
P.be2 = 1.0;%下层的等效刚度与基准刚度之比  P.be2 = 2*K2*(1-L)/L; 
P.mu  = 0.5;%下层质量对上层质量的比值
P.al1 = 0;  %上层准零刚度对应的等效线性刚度P.al1 = 2*K1*(1-L)/L; 
P.ga1 = 0.0;%上层准零刚度对应的三次方系数  P.ga1 = K1 / (U^2 * L^3);   
P.ga2 = 0.0;%下层准零刚度对应的三次方系数  P.ga2 = K2 / (U^2 * L^3); 
P.ze1 = 0.111803; %第二增阻尼比  需要根据质量比进行换算
P.lam   = 0.0;  %设置的第一层阻尼比为根据质量比进行计算，当电路断开时，模拟电磁分流阻尼比
P.kap_e = 0.0;  %电感系数
P.kap_c = 0.0;  %电容系数
P.sigma = 0.0;  %电阻系数
   
  
sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.05;          % 外激励幅值（小一些更线性）
FixedOmega = [];    % 扫频模式

fprintf('================================================\n');
fprintf('VERIFICATION 2: Linear Pure Mechanical Two-Stage (dB + Peak Picking)\n');
fprintf('================================================\n');
fprintf('Params: be1=%.3f, be2=%.3f, mu=%.3f, al1=%.3f, ze1=%.3f, lam=%.3f\n', ...
    P.be1, P.be2, P.mu, P.al1, P.ze1, P.lam);
fprintf('Circuit muted: kap_e=%.3f, kap_c=%.3f, sigma=%.3f\n', P.kap_e, P.kap_c, P.sigma);

%% -----------------------------
% 2) 理论线性无阻尼固有频率（经典 two-stage 2DOF）
% 说明：忽略阻尼与非线性，只用质量矩阵 M 和刚度矩阵 K
% -----------------------------
k12 = (P.be1 + P.al1);
k2  = P.be2;

M = diag([1, P.mu]);
K = [ k12,    -k12;
     -k12, k12 + k2];

wn = sqrt(eig(M\K));   % 无量纲固有频率
wn = sort(real(wn));

fprintf('\n--- Theoretical Undamped Natural Frequencies (2DOF) ---\n');
fprintf('  wn1 = %.6f\n', wn(1));
fprintf('  wn2 = %.6f\n', wn(2));
fprintf('------------------------------------------------------\n');

%% -----------------------------
% 3) 调用 FRF（HBM + 延拓扫频）
% -----------------------------
try
    x_results = FRF(sysP);
catch ME
    fprintf('Error in FRF: %s\n', ME.message);
    return;
end

if isempty(x_results)
    error('No results returned by FRF.');
end

%% -----------------------------
% 4) 数据处理 (转dB)
% -----------------------------
Omega_vec = x_results(16, :);

% 幅值 (Linear Scale) - 取 1次谐波幅值：sqrt(cos^2 + sin^2)
Amp_X1_Lin = sqrt(x_results(2,:).^2 + x_results(3,:).^2);
Amp_X2_Lin = sqrt(x_results(7,:).^2 + x_results(8,:).^2);

% 防止 log10(0) -> -Inf
eps_db = 1e-14;
Amp_X1_dB = 20 * log10(Amp_X1_Lin + eps_db);
Amp_X2_dB = 20 * log10(Amp_X2_Lin + eps_db);

%% -----------------------------
% 5) 自动寻峰（建议只抓前2个主峰）
% -----------------------------
opts.minProm_dB   = 1.0;      % 最小突起（dB），避免噪声点
opts.minDistOmega = 0.05;     % 峰间最小频率间距（按你扫频分辨率可改）
opts.maxNumPeaks  = 2;        % 2DOF 线性系统：只需要前2个主峰

[pks1, locs1] = find_peaks_custom(Amp_X1_dB, Omega_vec, opts);
[pks2, locs2] = find_peaks_custom(Amp_X2_dB, Omega_vec, opts);

fprintf('\n--- Resonance Peaks Report (Numerical FRF) ---\n');
fprintf('X1 (Upper Mass) Peaks:\n');
for i = 1:length(pks1)
    fprintf('  Peak %d: Omega = %.6f, Amp = %.2f dB\n', i, locs1(i), pks1(i));
end

fprintf('X2 (Lower Mass) Peaks:\n');
for i = 1:length(pks2)
    fprintf('  Peak %d: Omega = %.6f, Amp = %.2f dB\n', i, locs2(i), pks2(i));
end
fprintf('---------------------------------------------\n');

% 简单对照：把理论 wn 也打印到图上更“审稿友好”
%% -----------------------------
% 6) 绘图 (dB)
% -----------------------------
figure('Position',[100,100,1100,420], 'Color', 'w');

% X1
subplot(1,2,1);
plot(Omega_vec, Amp_X1_dB, 'b-', 'LineWidth', 2); hold on;
plot(locs1, pks1, 'ro', 'MarkerFaceColor', 'r'); 
xline(wn(1), '--k', 'wn1'); xline(wn(2), '--k', 'wn2');
xlabel('\Omega'); ylabel('|X_1| (dB)');
title('Linear Response: Upper Mass (m1)');
grid on; axis tight;

% X2
subplot(1,2,2);
plot(Omega_vec, Amp_X2_dB, 'r-', 'LineWidth', 2); hold on;
plot(locs2, pks2, 'bo', 'MarkerFaceColor', 'b');
xline(wn(1), '--k', 'wn1'); xline(wn(2), '--k', 'wn2');
xlabel('\Omega'); ylabel('|X_2| (dB)');
title('Linear Response: Lower Mass (m2)');
grid on; axis tight;

fprintf('Done! Peaks are marked with circles; wn1/wn2 are shown as dashed lines.\n');

%% ========================================================================
% 辅助函数：简单寻峰（不依赖工具箱）
% - 增加：最小突起 minProm_dB、峰间最小间距 minDistOmega、最多峰数 maxNumPeaks
% ========================================================================
function [pks, locs] = find_peaks_custom(y, x, opts)
    pks = []; locs = [];
    if length(y) < 3, return; end

    if ~isfield(opts, 'minProm_dB'),   opts.minProm_dB = 0; end
    if ~isfield(opts, 'minDistOmega'), opts.minDistOmega = 0; end
    if ~isfield(opts, 'maxNumPeaks'),  opts.maxNumPeaks = inf; end

    % 1) 先找所有局部极大
    cand_pk = [];
    cand_x  = [];
    for i = 2:length(y)-1
        if y(i) > y(i-1) && y(i) > y(i+1)
            cand_pk(end+1,1) = y(i); %#ok<AGROW>
            cand_x(end+1,1)  = x(i); %#ok<AGROW>
        end
    end
    if isempty(cand_pk), return; end

    % 2) 突起过滤：用“相邻谷值”近似（很轻量，不用工具箱）
    keep = false(size(cand_pk));
    for k = 1:length(cand_pk)
        % 左侧谷值
        il = find(x < cand_x(k), 1, 'last');
        if isempty(il), il = 1; end
        % 右侧谷值
        ir = find(x > cand_x(k), 1, 'first');
        if isempty(ir), ir = length(y); end

        left_min  = min(y(max(1,il-10):il));
        right_min = min(y(ir:min(length(y),ir+10)));

        prom = cand_pk(k) - max(left_min, right_min);
        if prom >= opts.minProm_dB
            keep(k) = true;
        end
    end
    cand_pk = cand_pk(keep);
    cand_x  = cand_x(keep);
    if isempty(cand_pk), return; end

    % 3) 按峰值从大到小排序
    [cand_pk, idx] = sort(cand_pk, 'descend');
    cand_x = cand_x(idx);

    % 4) 峰间最小距离筛选
    sel_pk = [];
    sel_x  = [];
    for k = 1:length(cand_pk)
        if isempty(sel_x)
            sel_pk(end+1,1) = cand_pk(k); %#ok<AGROW>
            sel_x(end+1,1)  = cand_x(k);  %#ok<AGROW>
        else
            if all(abs(cand_x(k) - sel_x) >= opts.minDistOmega)
                sel_pk(end+1,1) = cand_pk(k); %#ok<AGROW>
                sel_x(end+1,1)  = cand_x(k);  %#ok<AGROW>
            end
        end
        if length(sel_pk) >= opts.maxNumPeaks
            break;
        end
    end

    % 5) 最终按频率从小到大输出（更像“模态1、模态2”）
    [locs, idx2] = sort(sel_x, 'ascend');
    pks = sel_pk(idx2);
end
