function [x] = FRF(sysP)
%% FRF: HBM + 弧长延拓追踪幅频响应 (限制Omega范围版本)
% 输出 x: 16 x N
%  - 前15维: HBM 状态(主/副/电路的 0,1,3 次谐波系数)
%  - 第16维: Omega
%
% 关键改动：
% 1) 不再用 20000/30000 这种超大步数；
% 2) 增加 OmegaMax 终止条件，扫到目标频率立刻停止；
% 3) 结果自动裁剪为实际扫到的长度，避免扫到 200+

global Fw FixedOmega
FixedOmega = [];  % FRF 模式：nondim_temp2 把“第16维”当 Omega

fprintf('开始执行 FRF 扫频...\n');
tic;

%% ---------- 用户可调参数（建议你只改这里） ----------
Omega0   = 0.2;      % 起始频率
Omega1   = 0.21;     % 第二个点（用于初始切线）
OmegaMax = 10.0;     % !!! 关键：扫频上限（先10，后续可改30）
branch_max = 6000;   % 最大步数上限（配合 OmegaMax 使用即可）
print_every = 200;   % 打印频率

% 给一个微小初值避免雅可比奇异
y = zeros(15,1);
y(2) = 0.05;
% ------------------------------------------------------

%% ---------- 1) 构造两个初始点（Newton 固定 Omega） ----------
y_init0 = [y; Omega0];
x0_full = newton('nondim_temp2', y_init0, sysP);

y_init1 = [x0_full(1:15); Omega1];
x1_full = newton('nondim_temp2', y_init1, sysP);

x0 = x0_full(1:15);
x1 = x1_full(1:15);

%% ---------- 2) 预分配结果矩阵 ----------
% branch_follow2 通常会返回 (15+1) x N 或者 16 x N（看你实现）
% 这里我们自己做一个“截断器”：先跑到 branch_max，再裁剪到 OmegaMax
[x_all, ~] = branch_follow2('nondim_temp2', branch_max, Omega0, Omega1, x0, x1, sysP);

% 如果 branch_follow2 已经返回 16xN（最后一行为 Omega），直接处理；
% 如果返回 15xN（不含Omega），那说明你的 branch_follow2 内部没把Omega拼出来，需要你去 branch_follow2 里补。
if size(x_all,1) ~= 16
    error('FRF: branch_follow2 输出不是 16xN。请确认 branch_follow2 是否把 Omega 作为第16行返回。');
end

Omega_vec = x_all(end,:);

%% ---------- 3) 按 OmegaMax 截断 ----------
idx = find(Omega_vec <= OmegaMax);
if isempty(idx)
    warning('FRF: 没有扫到任何 Omega <= OmegaMax，检查起点或 newton 初值。');
    x = x_all;
else
    last = idx(end);
    x = x_all(:,1:last);
end

%% ---------- 4) 打印信息 ----------
last_omega = x(end,end);
fprintf('FRF 追踪完成。Omega 范围: %.3f -> %.3f\n', x(end,1), last_omega);
fprintf('FRF 点数: %d\n', size(x,2));
toc;

if last_omega < (OmegaMax*0.8)
    fprintf('[提示] 实际扫到的上限 %.2f 低于 OmegaMax=%.2f，可能在中途发生分支困难/收敛失败。\n', last_omega, OmegaMax);
end

end
