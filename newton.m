function [x, converged, R_norm] = newton(fun, x, sysP)
    %% 鲁棒版 Newton 求解器 (v2.0: 启用解析雅可比)
    % 修正点：直接调用 fun 的第二个输出 (Jac)，抛弃有限差分。
    
    max_iter = 50;      
    tol      = 1e-6;    % 适当放宽一点，1e-8 对 QZS 太苛刻
    
    converged = false;
    
    for k = 1:max_iter
        % --- 关键修改：同时获取残差 R 和 雅可比 J ---
        try
            if nargout(fun) >= 2
                [R, J_full] = feval(fun, x, sysP);
            else
                R = feval(fun, x, sysP);
                J_full = []; % 没提供解析解再说
            end
        catch
             % 兼容某些旧函数
             R = feval(fun, x, sysP);
             J_full = [];
        end
        
        n_eqs = length(R);
        n_vars = length(x);
        
        % 确定求解模式 (FRF vs Continuation)
        if n_vars == n_eqs
            active_idx = 1:n_vars; 
        elseif n_vars == n_eqs + 1
            active_idx = 1:n_eqs; % 最后一维是参数，固定
        else
            error('维度不匹配');
        end
        
        % 检查收敛
        R_norm = norm(R);
        if R_norm < tol, converged = true; return; end
        
        % --- 构造计算用的雅可比 J ---
        if ~isempty(J_full)
            % 【方案A】使用解析雅可比 (屠龙刀)
            J = J_full(:, active_idx);
        else
            % 【方案B】降级回有限差分 (指甲刀) - 仅作备用
            h_step = 1e-6;
            J = zeros(n_eqs, length(active_idx));
            for i = 1:length(active_idx)
                xt = x; xt(active_idx(i)) = xt(active_idx(i)) + h_step;
                Rt = feval(fun, xt, sysP);
                J(:, i) = (Rt - R) / h_step;
            end
        end
        
        % --- 求解修正量 (增加 Levenberg-Marquardt 正则化) ---
        % 当接近 QZS 时，J 可能奇异，加一点对角阻尼 lambda
        lambda = 1e-4 * norm(R); % 残差越大，阻尼越大，越像梯度下降
        J_damp = J' * J + lambda * eye(size(J,2));
        rhs    = -(J' * R);
        
        dx = J_damp \ rhs;
        
        % --- 简单的回溯线搜索 (防止发散) ---
        alpha = 1.0;
        for line_search = 1:5
            x_new = x;
            x_new(active_idx) = x_new(active_idx) + alpha * dx;
            R_new = feval(fun, x_new, sysP);
            if norm(R_new) < R_norm
                x = x_new;
                break;
            end
            alpha = alpha * 0.5;
        end
    end
end