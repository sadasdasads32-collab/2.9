function [x_branch, info] = L1(x_interest, sysP, step, branch_len, Omega_fixed)
%% L1  定频扫力
% 固定频率 Omega = FixedOmega (或 Omega_fixed)，以外力幅值 Fw 为延拓参数追踪周期解分支。
%
% 调用（强制对口）：
%   x_branch = L1(x_interest, sysP, step, branch_len)
%   x_branch = L1(x_interest, sysP, step, branch_len, Omega_fixed)
%
% x_interest:
%   - 15x1 : 15个HBM系数，初始力 F0 从 global Fw 读取
%   - 16x1 : [15个HBM系数; F0]  ★第16个是 F0（不是Omega）
%
% 输出：
%   x_branch : 16xN, 前15行HBM系数，最后1行是F（延拓参数）
%   info     : 结构体，记录Omega、F0、步长等信息
%
% 依赖：
%   newton.m, branch_follow2.m, branch_aux2.m, nondim_temp2.m
%
% 与 nondim_temp2 的对口约定：
%   - 当 global FixedOmega 非空：nondim_temp2 将 y(16) 解释为 Fw（扫力）
%   - 当 global FixedOmega 为空：nondim_temp2 将 y(16) 解释为 Omega（扫频）

    global Fw FixedOmega

    %% -----------------------------
    % 0) 默认值与健壮性
    % -----------------------------
    if nargin < 3 || isempty(step),       step = 5e-4;    end
    if nargin < 4 || isempty(branch_len), branch_len = 4000; end

    if ~isscalar(step) || ~isfinite(step) || step == 0
        error('L1: step must be a finite nonzero scalar.');
    end
    if ~isscalar(branch_len) || ~isfinite(branch_len) || branch_len < 50
        error('L1: branch_len must be a finite scalar >= 50.');
    end
    branch_len = round(branch_len);

    % 防止步长太小卡死
    step_min = 1e-8;
    if abs(step) < step_min
        step = sign(step) * step_min;
    end

    %% -----------------------------
    % 1) 解析固定频率 Omega（严格）并强制进入扫力模式
    % -----------------------------
    FixedOmega_backup = FixedOmega;
    Fw_backup         = Fw;

    if nargin >= 5 && ~isempty(Omega_fixed)
        if ~isscalar(Omega_fixed) || ~isfinite(Omega_fixed) || Omega_fixed <= 0
            error('L1: Omega_fixed must be a positive finite scalar.');
        end
        FixedOmega = Omega_fixed;
    else
        if isempty(FixedOmega) || ~isfinite(FixedOmega) || FixedOmega <= 0
            error(['L1: FixedOmega is empty/invalid. ' ...
                   'Set global FixedOmega before calling, or pass Omega_fixed as 5th input.']);
        end
        % 使用已有 FixedOmega
    end

    Omega_use = FixedOmega;

    %% -----------------------------
    % 2) 解析初值：HBM系数 + 初始力 F0
    % -----------------------------
    x_interest = x_interest(:);

    if numel(x_interest) < 15
        error('L1: x_interest must have at least 15 elements. Got %d.', numel(x_interest));
    end

    x0_coeff = x_interest(1:15);

    if numel(x_interest) >= 16
        % ★约定：第16个是初始力 F0（不是Omega）
        F0 = x_interest(16);
    else
        if isempty(Fw) || ~isfinite(Fw)
            error('L1: global Fw is empty/invalid, cannot infer initial force F0.');
        end
        F0 = Fw;
    end

    if ~isfinite(F0)
        error('L1: initial force F0 is NaN/Inf.');
    end

    % ★同步全局Fw，避免 nondim_temp2 里仍使用 global Fw 造成不一致
    Fw = F0;

    %% -----------------------------
    % 3) 两点初始化 + Newton校正（强烈推荐）
    % -----------------------------
    y0_guess = [x0_coeff; F0];
    y1_guess = [x0_coeff; F0 + step];

    ok0 = false; ok1 = false;

    %% -----------------------------
    % 3) 两点初始化 + Newton校正（强烈推荐）
    % -----------------------------
    y0_guess = [x0_coeff; F0];
    
    % 获取 newton 的第二个输出作为成功标志
    [y0, ok0] = newton('nondim_temp2', y0_guess, sysP);
    if ~ok0
        warning('L1: 第1个点 Newton 未收敛，可能偏离解流形！ F0=%.6g', F0);
        y0 = y0_guess; % 强行沿用猜测值
    end
    
    % 用 y0(1:15) 作为 y1 的系数种子，更稳
    y1_guess = [y0(1:15); F0 + step];
    
    [y1, ok1] = newton('nondim_temp2', y1_guess, sysP);
    if ~ok1
        warning('L1: 第2个点 Newton 未收敛！强制继续可能会导致延拓飞点。 F1=%.6g', F0+step);
        y1 = y1_guess;
    end
    
    mu0 = y0(end);   % 扫力参数：F0
    mu1 = y1(end);   % 扫力参数：F0+step
    fprintf('L1(force): Omega=%.6f | F0=%.6g -> F1=%.6g | step=%.2e | len=%d | Newton[%d,%d]\n', ...
            Omega_use, mu0, mu1, step, branch_len, ok0, ok1);

    % 用y0(1:15)作为y1的系数种子，更稳
    y1_guess = [y0(1:15); F0 + step];

    try
        % 同步全局Fw为第二点力（可显著提升收敛一致性）
        Fw = F0 + step;
        y1 = newton('nondim_temp2', y1_guess, sysP);
        ok1 = true;
    catch ME
        warning('L1: Newton failed at F1=%.6g | %s', F0+step, ME.message);
        y1 = y1_guess;
    end

    mu0 = y0(end);   % 扫力参数：F0
    mu1 = y1(end);   % 扫力参数：F0+step

    fprintf('L1(force): Omega=%.6f | F0=%.6g -> F1=%.6g | step=%.2e | len=%d | Newton[%d,%d]\n', ...
            Omega_use, mu0, mu1, step, branch_len, ok0, ok1);

    %% -----------------------------
    % 4) 弧长延拓：参数是力 F
    % -----------------------------
    % ★延拓期间 global Fw 不再固定，它实际由 y(16)=F 控制；
    % 但为了兼容某些实现，我们把 global Fw 设置为起点力，通常足够稳定。
    Fw = mu0;

    [x_branch, conv] = branch_follow2('nondim_temp2', branch_len, mu0, mu1, y0(1:15), y1(1:15), sysP);

    if ~isempty(conv)
        % 这里不强依赖 conv 的语义（不同工程可能不同）
    end

    fprintf('L1(force): done | N=%d | F_end=%.6g\n', size(x_branch,2), x_branch(end,end));

    %% -----------------------------
    % 5) 输出 info
    % -----------------------------
    info = struct();
    info.mode = 'force_continuation';
    info.Omega_fixed = Omega_use;
    info.F0 = mu0;
    info.F1 = mu1;
    info.step = step;
    info.branch_len = branch_len;
    info.newton_ok = [ok0, ok1];

    %% -----------------------------
    % 6) 恢复全局变量（避免污染其他脚本）
    % -----------------------------
    FixedOmega = FixedOmega_backup;
    Fw         = Fw_backup;
end
