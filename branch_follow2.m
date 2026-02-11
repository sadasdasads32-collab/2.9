function [x, p] = branch_follow2(fname, nsteps, mu0, mu1, x0, x1, sysP)
% 弧长连续法路径追踪 (审稿人最终修正版：物理范围限制 2.0)
% 功能：自适应步长 + 切向预测 + 物理范围锁死

    global tracking_file_name xc arc 
    % 注意：必须保留 xc 和 arc 为全局变量，因为 branch_aux2 需要用到它们

    tracking_file_name = fname;

    % --- 1. 数据预处理 ---
    x0 = x0(:);
    x1 = x1(:);
    
    % 拼装成完整状态向量 [State; Parameter]
    if numel(x0) == 15, x0 = [x0; mu0]; end
    if numel(x1) == 15, xc = [x1; mu1]; end
    
    % 如果输入本身就是 16 维，则直接使用
    if numel(x0) == 16, x0 = x0; end 
    if numel(xc) == 16, xc = xc; end 

    % 初始化结果存储
    x = [x0, xc];
    
    % --- 2. 初始步长设置 ---
    % 计算初始切向向量
    tangent = xc - x0;
    current_arc = norm(tangent); 
    if current_arc < 1e-8
        current_arc = 1e-3; % 防止起点重合
        tangent = [zeros(15,1); 1e-3]; 
    end
    tangent = tangent / norm(tangent); % 归一化切向
    
    % 同步全局变量 arc (供 newton -> branch_aux2 使用)
    arc = current_arc; 
    
    % 步长控制参数
    arc_min = 1e-5;
    arc_max = 0.1;    % 最大步长
    opt_iter = 4;     % 最佳迭代次数

    k = 1;
    p = 'y'; 
    
    fprintf('Branch Follow (Final): Start. Arc=%.1e, Limit=2.0\n', arc);

    while k < nsteps
        
        % --- 3. 预测 (Predictor): 基于切向 ---
        xg = xc + tangent * arc; 
        
        % --- 4. 校正 (Corrector): 牛顿迭代 ---
        [xx, ok, Rn, iter_count] = try_solve(fname, xg, sysP);
        
        % --- 5. 步长自适应策略 (Adaptive Logic) ---
        if ok
            % === 收敛成功 ===
            k = k + 1;
            
            % 更新切向
            new_tangent = xx - xc;
            dist = norm(new_tangent);
            if dist > 1e-12
                tangent = new_tangent / dist;
            end
            
            % 存数
            x0 = xc;
            xc = xx;
            x = [x, xx];
            
            % 进度打印 (每50步)
            if mod(k, 50) == 0
                fprintf('   Step %4d: Force=%.4f | Arc=%.1e | Iter=%d\n', ...
                        k, xc(end), arc, iter_count);
            end
            
            % >> 步长调整 <<
            if iter_count < opt_iter
                arc = min(arc * 1.5, arc_max); % 加速
            elseif iter_count > opt_iter + 2
                arc = max(arc * 0.7, arc_min); % 减速
            end
            
            % >> 物理范围检查 (User Requested Fix) <<
            % 将上限锁死在 2.0，超过即停
            if xc(end) < -0.05 || xc(end) > 2.0 
                fprintf('   End: Parameter out of range (reached limit %.2f).\n', xc(end));
                break;
            end
            
        else
            % === 收敛失败 (Cut Step) ===
            fprintf('   Step %d failed (Res=%.2e). Retrying with smaller arc...\n', k, Rn);
            arc = arc * 0.5; % 步长减半重试
            
            if arc < arc_min
                fprintf('   Error: Arc length too small (%.1e). Stopping.\n', arc);
                p = 'n';
                break;
            end
            continue; 
        end
    end
    
    fprintf('Branch Follow: Finished %d steps.\n', k);
end

% --- 辅助子函数 ---
function [xx, ok, Rn, iter] = try_solve(fname, xg, sysP)
    [xx, ok, Rn] = newton('branch_aux2', xg, sysP);
    
    % 估算迭代次数用于步长控制
    if Rn < 1e-7
        iter = 3; 
    elseif Rn < 1e-4
        iter = 5; 
    else
        iter = 10; 
    end
    
    if Rn > 1e-3, ok = 0; end
end