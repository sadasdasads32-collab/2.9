function [x, converged, R_norm] = newton(fun, x, sysP)
    %% 鲁棒版 Newton 求解器 (通用型)
    % 自动处理 FRF (16输入/15输出) 和 分岔追踪 (31输入/31输出) 两种模式
    % 解决了 "需要雅可比矩阵" 和 "矩阵奇异" 等问题
    
    % 说明：
    % - 第 2 个输出 converged 为逻辑量，用于 branch_follow2 的收敛判断。
    % - 第 3 个输出 R_norm 为最终残差范数，便于 residual 检查/作图。
    %
    % 参数（可按需要再收紧/放松）
    max_iter = 80;      % 最大迭代次数
    tol      = 1e-8;    % 收敛容差（建议 <=1e-6，做论文级 residual 用 1e-8 更稳）
    h_step   = 1e-5;    % 数值微分步长（太大影响精度，太小放大噪声）

    converged = false;
    
    for k = 1:max_iter
        % 1. 计算当前残差
        z = feval(fun, x, sysP);
        R = z; % 假定输出 z 就是残差向量
        
        n_eqs = length(R);
        n_vars = length(x);
        
        % 2. 智能判断求解模式
        if n_vars == n_eqs
            % [模式 A] 方阵系统 (如 L1 分岔追踪: 31x31)
            % 所有变量都参与迭代
            active_indices = 1:n_vars;
            
        elseif n_vars == n_eqs + 1
            % [模式 B] 参数扫描 (如 FRF 扫频: 16输入 -> 15方程)
            % 最后一个变量是参数(频率)，固定不动，只迭代前 n_eqs 个变量
            active_indices = 1:n_eqs;
            
        else
            error('Newton: 维度异常！输入变量数(%d)与方程数(%d)不匹配。', n_vars, n_eqs);
        end
        
        % 3. 检查收敛性
        R_norm = norm(R);
        if R_norm < tol
            converged = true;
            return; % 收敛，返回结果
        end
        
        % 4. 数值计算雅可比矩阵 J (只对有效变量求导)
        J = zeros(n_eqs, length(active_indices));
        
        for i = 1:length(active_indices)
            idx = active_indices(i);
            
            % 施加微扰
            x_temp = x;
            x_temp(idx) = x_temp(idx) + h_step;
            
            % 计算扰动后的残差
            z_temp = feval(fun, x_temp, sysP);
            
            % 差分求导
            J(:, i) = (z_temp - R) / h_step;
        end
        
        % 5. 求解修正量 dx
        % 使用左除 (\) 代替 inv()，即使矩阵接近奇异也能求出最优解
        if cond(J) > 1e12
            % 矩阵极度病态时，使用阻尼最小二乘法 (Levenberg-Marquardt 策略)
            dx = -(J' * J + 0.01 * eye(size(J))) \ (J' * R);
        else
            % 正常情况使用高斯消元
            dx = -J \ R;
        end
        
        % 6. 阻尼/线搜索更新（避免折叠附近残差不降）
        alpha = 1.0;
        x_trial = x;
        x_trial(active_indices) = x_trial(active_indices) + alpha * dx;
        R_trial = feval(fun, x_trial, sysP);
        Rn_trial = norm(R_trial);

        % 如果残差没有下降，则逐步减小步长
        ls_iter = 0;
        while (Rn_trial > R_norm) && (alpha > 1e-4) && (ls_iter < 12)
            alpha = alpha * 0.5;
            x_trial = x;
            x_trial(active_indices) = x_trial(active_indices) + alpha * dx;
            R_trial = feval(fun, x_trial, sysP);
            Rn_trial = norm(R_trial);
            ls_iter = ls_iter + 1;
        end

        x = x_trial;
        R_norm = Rn_trial;
    end
    
    % 如果循环结束还没收敛，返回 converged=false，同时输出最终残差。
end