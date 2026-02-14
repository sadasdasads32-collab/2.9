function [x, p] = branch_follow2(fname, nsteps, mu0, mu1, x0, x1, sysP)
% branch_follow2: 弧长连续法路径追踪（强稳增强版）
% 改进点：
% 1) 弧长范数加权：参数维度权重 wp，避免状态尺度淹没参数
% 2) 收敛判据更合理：优先相信 newton 的 ok，Rn 只做辅助
% 3) 多次失败后自动“切向重置”为偏参数方向，提高起步成功率
% 4) retry+自适应步长，支持 FRF/L1

    global tracking_file_name xc arc
    global ParamMin ParamMax

    tracking_file_name = fname;

    %% ---------- 1) 输入预处理 ----------
    x0 = x0(:); x1 = x1(:);
    if numel(x0)==15, x0=[x0; mu0]; end
    if numel(x1)==15, x1=[x1; mu1]; end
    if numel(x0)~=16 || numel(x1)~=16
        error('branch_follow2: x0/x1 must be 15x1 or 16x1.');
    end
    xc = x1;

    %% ---------- 2) 参数范围 ----------
    if isempty(ParamMin) || ~isfinite(ParamMin), ParamMin = min(mu0,mu1) - 0.05; end
    if isempty(ParamMax) || ~isfinite(ParamMax), ParamMax = max(mu0,mu1) + 0.05; end

    dir = sign(mu1-mu0); if dir==0, dir=1; end

    %% ---------- 3) 弧长/切向初始化（加权） ----------
    wp = 10;  % ★参数维度权重（关键参数！） 5~20 之间可调
    tangent = (xc - x0);

    arc0 = weighted_norm(tangent, wp);
    if arc0 < 1e-12
        tangent = [zeros(15,1); dir];
        arc0 = 1e-3;
    end
    tangent = tangent / max(weighted_norm(tangent, wp), 1e-12);

    arc = arc0;

    %% ---------- 4) 步长控制 ----------
    arc_min    = 1e-6;
    arc_max    = 5e-2;
    arc_grow   = 1.25;
    arc_shrink = 0.5;

    retry_max  = 15;
    fail_reset_after = 6;  % 连续失败这么多次，就重置切向

    %% ---------- 5) 预分配 ----------
    x = zeros(16, nsteps+2);
    x(:,1)=x0; x(:,2)=xc;
    col=2; k=1; p='y';

    fprintf('Branch Follow (Weighted): Start. Arc=%.2e, wp=%g, Range=[%.6f, %.6f], dir=%+d\n',...
        arc, wp, ParamMin, ParamMax, dir);

    consec_fail = 0;

    %% ---------- 6) 主循环 ----------
    while k < nsteps

        success = false;
        arc_try = arc;

        for rtry = 1:retry_max
            % predictor（注意：对参数范围先检查）
            xg = xc + tangent * arc_try;

            if xg(end) < ParamMin || xg(end) > ParamMax
                arc_try = arc_try * arc_shrink;
                if arc_try < arc_min, break; end
                continue;
            end

            [xx, ok, Rn, iter_est] = try_solve(xg, sysP);

            if ok
                success = true;
                break;
            else
                arc_try = arc_try * arc_shrink;
                if arc_try < arc_min, break; end
            end
        end

        if ~success
            consec_fail = consec_fail + 1;

            % 连续失败太多：重置切向为“偏参数方向”，救起步
            if consec_fail >= fail_reset_after
                tangent = [zeros(15,1); dir];
                tangent = tangent / weighted_norm(tangent, wp);
                consec_fail = 0;
                arc = max(arc_min*10, arc*0.2);
                fprintf('   [Reset] tangent -> param direction, arc -> %.2e\n', arc);
                continue;
            end

            fprintf('   Fail: step=%d, arc_try=%.2e < arc_min or no convergence. Stop.\n', k, arc_try);
            p='n'; x=x(:,1:col); return;
        end

        % accept
        consec_fail = 0;
        k = k + 1; col = col + 1;
        x(:,col) = xx;

        % 更新切向（加权归一）
        new_tan = xx - xc;
        dn = weighted_norm(new_tan, wp);
        if dn > 1e-14
            tangent = new_tan / dn;
        end

        % 更新点
        x0 = xc; xc = xx;

        % 用本次成功的 arc_try 作为基准更新 arc
        arc = arc_try;

        % 步长自适应（iter_est 只是粗估）
        if iter_est <= 5
            arc = min(arc*arc_grow, arc_max);
        elseif iter_est >= 10
            arc = max(arc*0.8, arc_min);
        end

        if mod(k,50)==0 || k<=5
            fprintf('   Step %4d: Param=%.6f | arc=%.2e | R=%.2e | it~%d\n',...
                k, xc(end), arc_try, Rn, iter_est);
        end

        if xc(end) < ParamMin || xc(end) > ParamMax
            fprintf('   End: reached param bound %.6f (Range=[%.6f, %.6f])\n',...
                xc(end), ParamMin, ParamMax);
            break;
        end
    end

    x = x(:,1:col);
    fprintf('Branch Follow (Weighted): Finished. Steps=%d, LastParam=%.6f\n', k, x(end,end));

end

%% ====== 加权范数：状态+参数权重 ======
function n = weighted_norm(v, wp)
    ds = v(1:15);
    dp = v(16);
    n = sqrt(sum(ds.^2) + (wp*dp)^2);
end

%% ====== Newton 调用封装（只依赖 newton 返回 [xx, ok, Rn]）=====
function [xx, ok, Rn, iter_est] = try_solve(xg, sysP)
    [xx, ok, Rn] = newton('branch_aux2', xg, sysP);

    if ~isfinite(Rn)
        ok = 0; Rn = inf; iter_est = 99; return;
    end

    % iter_est 粗估（用于步长调节，不作为生死判据）
    if Rn < 1e-10
        iter_est = 3;
    elseif Rn < 1e-8
        iter_est = 4;
    elseif Rn < 1e-6
        iter_est = 6;
    elseif Rn < 1e-4
        iter_est = 9;
    else
        iter_est = 12;
    end

    % ★关键：不再用 “Rn>5e-4 就直接判死刑”
    % 只要 newton 给 ok=1，就先接受，让延拓推进（特别是起步阶段）
    % 你若担心垃圾解，可在这里加：if Rn>1e-2, ok=0; end
    if Rn > 1e-2
        ok = 0;
    end
end
