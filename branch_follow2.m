function [x,p]=branch_follow2(fname,nsteps,mu0,mu1,x0,x1,sysP)
% 弧长连续法路径追踪 (自动化修正版)
% 移除了所有手动交互 (input)，防止死循环
% 移除了 mu > 1 的硬性限制

global mu tracking_file_name xc arc

% 指定用于构建弧长约束方程的文件名 (通常是 branch_aux2)
tracking_file_name=fname;

% 扩展状态向量：[变量 x; 参数 mu]
% 兼容两种输入：
%  - 仅给 15 个系数（此时这里追加 mu）
%  - 已给 16 维（15系数 + mu），此时忽略 mu0/mu1 形参
x0 = x0(:);
x1 = x1(:);
if numel(x0) == 16
    mu0 = x0(end);
    x0 = x0; % already includes mu
else
    x0 = [x0; mu0];
end

if numel(x1) == 16
    mu1 = x1(end);
    xc = x1;
else
    xc = [x1; mu1];
end

% 初始化存储矩阵
x=[x0,xc];

% 根据前两点确定弧长步长 (固定)
arc=norm(x0-xc);

k=1;
c=1;
p='y'; % 默认为继续

fprintf('Branch Follow: 开始追踪... 目标步数: %d, 弧长: %.4e\n', nsteps, arc);

while (k<nsteps)*c
   % 1. 预测 (Predictor): 线性外推
   xg=2*xc-x0;  
   
   % 2. 校正 (Corrector): 牛顿迭代
   % 调用 newton 解 branch_aux2 (即 原方程 + 弧长约束)
   % 注意：newton 的第 2 个输出是逻辑收敛标志，第 3 个输出是残差范数
   [xx, ok, Rn] = newton('branch_aux2', xg, sysP);
   
   % 2.1 收敛性与残差门限判断
   % 论文级 residual 建议 <= 1e-4；如果你的 nondim_temp2 标度很大可适当放宽
   res_tol = 1e-4;
   if ok && (Rn < res_tol)
      % --- 收敛成功 ---
      k=k+1;
      
      % 更新点位
      x0=xc; % 旧点
      xc=xx; % 新点
      x=[x,xx]; % 存入结果
      
      % 实时保存结果 (可选：如果嫌慢可以注释掉下面这行)
      % save result x; 
      
      % --- 进度打印 ---
      % 每 100 步打印一次，防止刷屏
      if mod(k, 100) == 0
          fprintf('   Step %d/%d done. Parameter(mu) = %.4f\n', k, nsteps, xc(end));
      end
      % --- 新增：物理范围限制（这里 mu 表示外力幅值 F） ---
      current_mu = xc(end);
      mu_min = 0;      % 力幅值下限
      mu_max = 300;    % 力幅值上限（按需要调整）
      if current_mu > mu_max || current_mu < mu_min
          fprintf('Branch Follow: 力 F = %.2f 超出物理范围 [%.2f, %.2f]，停止追踪。\n', current_mu, mu_min, mu_max);
          break;
      end
      % --- 移除原有的 limit point 暂停逻辑 ---
      % 原代码检测到 abs(x0(end)-xc(end)) 很小时会暂停，现已移除，改为自动继续。
      
      % --- 移除原有的 mu > 1 暂停逻辑 ---
      % 你的频率要去到 50，不能在这里限制 > 1。
      
   else
       % --- 不收敛 / 残差过大：尝试自动缩小弧长再校正 ---
       shrink_ok = false;
       arc0 = arc;
       for shrink = 1:6
           arc = arc0 * 0.5^shrink; % 缩弧长
           [xx2, ok2, Rn2] = newton('branch_aux2', xg, sysP);
           if ok2 && (Rn2 < res_tol)
               xx = xx2; ok = ok2; Rn = Rn2;
               shrink_ok = true;
               break;
           end
       end

       if shrink_ok
           k = k + 1;
           x0 = xc;
           xc = xx;
           x = [x, xx];
           if mod(k, 100) == 0
               fprintf('   Step %d/%d done. Parameter(mu) = %.4f (res=%.2e, arc=%.2e)\n', k, nsteps, xc(end), Rn, arc);
           end
       else
           fprintf('Branch Follow: 在第 %d 步校正失败 (ok=%d, res=%.2e)，追踪终止。\n', k, ok, Rn);
           p = 'n';
           break;
       end
   end
end

fprintf('Branch Follow: 追踪完成。共计算 %d 步。\n', k);

end