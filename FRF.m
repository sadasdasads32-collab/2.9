function [x] = FRF(sysP)
%% FRF: HBM + 弧长延拓追踪幅频响应（强稳版：自适应起步 + 自动重启）
% 输出 x: 16 x N
%  - 1:15 为 HBM 状态（0/1/3 次谐波系数）
%  - 16 为 Omega
%
% 依赖：newton.m, branch_follow2.m, nondim_temp2.m, branch_aux2.m
%
% 关键思想：
% 1) 用 very small dOmega 构造第2个点，保证初始弧长不要过大
% 2) 如果 Newton 在 Omega1 失败，自动缩小 dOmega 重试
% 3) 如果 branch_follow2 早停，自动用最后两点重新起步继续扫（拼接）

global Fw FixedOmega ParamMin ParamMax

FixedOmega = [];  % FRF 模式：nondim_temp2 把第16维当 Omega

fprintf('开始执行 FRF 扫频 (robust)...\n');
tic;

%% ---------- 用户可调参数 ----------
Omega0   = 0.2;        % 起点
OmegaMax = 10.0;       % 上限
branch_max_total = 8000; % 总点数预算（允许自动重启拼接）
chunk_steps = 1200;    % 每次延拓的最大步数（分段更稳）
% 起步 dOmega（会自适应缩小）
dOmega_try = 1e-2;     % 初始尝试步长（不行会砍半）
dOmega_min = 1e-5;     % 最小起步步长（再小就不玩了）
% 初值扰动（越小越线性、越稳）
init_eps = 1e-4;       % 你原来 0.01 太大了
% 允许的最小有效点数（少于这个认为“失败需要重启”）
min_pts_ok = 20;
% -----------------------------------

ParamMin = Omega0 - 0.05;
ParamMax = OmegaMax;

%% ---------- 0) 第一个点：Newton @ Omega0 ----------
y0 = zeros(15,1);
y0(2) = init_eps;  % 轻微扰动避免奇异（也可全0）
y_init0 = [y0; Omega0];

[x0_full, ok0] = newton_safe('nondim_temp2', y_init0, sysP);
if ~ok0
    error('FRF: Newton 在 Omega0=%.6f 失败。请减小非线性/耦合或换更小 init_eps。', Omega0);
end
x0 = x0_full(1:15);

%% ---------- 1) 自适应构造第二点：Newton @ Omega1 ----------
dOm = dOmega_try;
ok1 = false;
x1_full = [];

while dOm >= dOmega_min
    Omega1 = Omega0 + dOm;
    y_init1 = [x0; Omega1];              % 用上一个解当初值最稳
    [x1_full, ok1] = newton_safe('nondim_temp2', y_init1, sysP);
    if ok1
        break;
    end
    dOm = dOm * 0.5;
end

if ~ok1
    error('FRF: 无法构造第二点（dOmega 已降到 %.1e 仍失败）。建议降低 lam/kap_e/非线性或先用更小 Fw。', dOmega_min);
end

x1 = x1_full(1:15);
Omega1 = Omega0 + dOm;

fprintf('起步成功：Omega0=%.6f, Omega1=%.6f (dOm=%.2e)\n', Omega0, Omega1, dOm);

%% ---------- 2) 分段延拓 + 自动重启拼接 ----------
x_all = [];
cur_x0 = x0;
cur_x1 = x1;
cur_Om0 = Omega0;
cur_Om1 = Omega1;

remain_budget = branch_max_total;

while remain_budget > 0
    nsteps = min(chunk_steps, remain_budget);

    fprintf('\n[FRF] 延拓分段：start=%.6f -> %.6f, nsteps=%d\n', cur_Om0, cur_Om1, nsteps);

    try
        [x_seg, conv] = branch_follow2('nondim_temp2', nsteps, cur_Om0, cur_Om1, cur_x0, cur_x1, sysP);
    catch ME
        fprintf('FRF: branch_follow2 崩溃：%s\n', ME.message);
        break;
    end

    % 基本检查
    if isempty(x_seg) || size(x_seg,1) ~= 16 || size(x_seg,2) < 2
        fprintf('FRF: 本段输出无效，停止。\n');
        break;
    end

    % 截掉 NaN/Inf
    Om = x_seg(16,:);
    good = isfinite(Om);
    if ~all(good)
        last_good = find(good, 1, 'last');
        if isempty(last_good) || last_good < 2
            fprintf('FRF: 本段一开始就 NaN/Inf，停止。\n');
            break;
        end
        x_seg = x_seg(:,1:last_good);
        Om = x_seg(16,:);
        fprintf('[诊断] 本段出现 NaN/Inf，截断到 Omega=%.6f\n', Om(end));
    end

    % 拼接（避免重复第1点）
    if isempty(x_all)
        x_all = x_seg;
    else
        x_all = [x_all, x_seg(:,2:end)]; %#ok<AGROW>
    end

    % 是否到达上限
    if x_all(16,end) >= OmegaMax
        fprintf('已到达 OmegaMax=%.3f，结束。\n', OmegaMax);
        break;
    end

    % 如果本段点数太少，说明延拓卡死：尝试“更小的起步间隔”重启
    if size(x_seg,2) < min_pts_ok
        fprintf('[FRF] 本段点数=%d 太少，尝试缩小起步间隔重启...\n', size(x_seg,2));
        % 用最后一个解作为新起点，并尝试很小 dOm 构造第二点
        cur_Om0 = x_all(16,end);
        cur_x0  = x_all(1:15,end);

        dOm2 = 5e-3;
        ok2 = false;
        while dOm2 >= dOmega_min
            cur_Om1 = cur_Om0 + dOm2;
            [tmp, ok2] = newton_safe('nondim_temp2', [cur_x0; cur_Om1], sysP);
            if ok2, break; end
            dOm2 = dOm2 * 0.5;
        end
        if ~ok2
            fprintf('重启失败：无法在 Omega=%.6f 附近构造第二点，停止。\n', cur_Om0);
            break;
        end
        cur_x1 = tmp(1:15);
        continue;
    end

    % 正常情况下，用最后两点作为下一段的起步（保证切向合理）
    cur_Om0 = x_all(16,end-1);
    cur_x0  = x_all(1:15,end-1);

    cur_Om1 = x_all(16,end);
    cur_x1  = x_all(1:15,end);

    remain_budget = branch_max_total - size(x_all,2);
end

%% ---------- 3) 最终按 OmegaMax 截断 ----------
Om_all = x_all(16,:);
idx = find(Om_all <= OmegaMax);
if isempty(idx)
    x = x_all;
else
    x = x_all(:,1:idx(end));
end

fprintf('\nFRF 完成：Omega %.6f -> %.6f, 点数=%d\n', x(16,1), x(16,end), size(x,2));
toc;

end

%% ========== 安全 Newton 封装（兼容你 newton 返回格式）==========
function [x, ok] = newton_safe(funname, x0, sysP)
    ok = false; x = [];
    try
        out = newton(funname, x0, sysP);
        % 兼容：有的 newton 返回 [x, ok]，有的只返回 x
        if iscell(out)
            x = out{1}; ok = out{2};
        else
            x = out;
            ok = true;
        end
    catch
        ok = false;
        x = [];
    end
end
