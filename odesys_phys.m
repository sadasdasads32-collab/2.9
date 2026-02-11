function dydt = odesys_phys(t, y, sysP, Omega, F_ampl)
    % odesys_phys: 时域验证用的标准 ODE 方程
    % 对应 nondim_temp2.m 的物理设置
    %
    % 状态向量 y (6x1): [x1; x1_dot; x2; x2_dot; q; q_dot]
    
    % --- 1. 解包参数 ---
    be1 = sysP(1); be2 = sysP(2); mu_mass = sysP(3);
    al1 = sysP(4); ga1 = sysP(5); ze1 = sysP(6);
    lam_phys = sysP(7); 
    kap_e = sysP(8); kap_c = sysP(9); sigma = sysP(10); ga2 = sysP(11);
    
    theta = sqrt(max(0, lam_phys)); % 对称耦合系数
    
    % --- 2. 状态变量 ---
    x1 = y(1); v1 = y(2);
    x2 = y(3); v2 = y(4);
    q  = y(5); i  = y(6); % i = q_dot
    
    % --- 3. 中间变量 ---
    x12 = x1 - x2;
    v12 = v1 - v2;
    
    % 恢复力 (Restoring Forces)
    % 上层 QZS: F_upper = (be1+al1)*x12 + ga1*x12^3
    F_upper = (be1 + al1) * x12 + ga1 * (x12^3);
    
    % 下层支撑: F_lower = be2*x2 + ga2*x2^3
    F_lower = be2 * x2 + ga2 * (x2^3);
    
    % 下层阻尼: F_damp2 = 2*mu*ze1 * v2 (注意归一化系数)
    F_damp2 = (2 * mu_mass * ze1) * v2; 
    
    % 电磁力 (EM Force)
    % 注意：根据 nondim_temp2，R1 中是 +force_em，且 force_em = theta*W*q' = theta*i
    % R2 中是 -force_em
    F_em = theta * i; 
    
    % --- 4. 动力学方程 ---
    % m1*x1'' + F_upper + F_em = F*cos(Wt)  (m1=1)
    % => x1'' = F*cos(Wt) - F_upper - F_em
    dx1 = v1;
    dv1 = F_ampl * cos(Omega * t) - F_upper - F_em;
    
    % m2*x2'' + F_damp2 + F_lower - F_upper - (-F_em) = 0 (m2=mu)
    % => mu*x2'' = F_upper - F_lower - F_damp2 - F_em
    % (注意：R2 中的 -force_em 移到等式右边变 +，但原方程 R2 = Inertia + ... - force_em = 0
    %  => Inertia = force_em ... 
    %  Wait, let's verify signs from nondim_temp2 strictly:
    %  R1 = Inertia1 + Stiffness1 + theta*i - Force = 0 => x1'' = Force - Stiffness1 - theta*i (Correct)
    %  R2 = Inertia2 + Damping2 + Stiffness2 - Stiffness1 - theta*i = 0
    %  => mu*x2'' = Stiffness1 + theta*i - Stiffness2 - Damping2 (Correct)
    dx2 = v2;
    dv2 = (F_upper + F_em - F_lower - F_damp2) / mu_mass;
    
    % Circuit: L*q'' + R*q' + 1/C*q - theta*(x1'-x2') = 0
    % => kap_e*i' + sigma*i + kap_c*q - theta*v12 = 0
    dq = i;
    di = (theta * v12 - sigma * i - kap_c * q) / kap_e;
    
    dydt = [dx1; dv1; dx2; dv2; dq; di];
end