function [z, Jac] = nondim_temp2(y, sysP)
%% nondim_temp2：T1（标准 Two-stage）物理一致性最终修正版
%
% 【修正历史】
%   1. 实施对称耦合：系数统一为 theta = sqrt(lam)
%   2. 修正电路源项符号：R3 = LHS - (+Source)
%   3. 【关键】修正机械阻尼符号：
%      - R1 中为 +force_em (正阻尼)
%      - R2 中为 -force_em (反作用力)
%      - 雅可比 J13, J23 同步修正
%
% 【输入】
%   y   : 16x1 [x1(5); x2(5); q(5); Omega_or_Fw]
%   sysP: 11x1 参数
%
% 【输出】
%   z   : 15x1 残差
%   Jac : 15x15 雅可比

    global Fw FixedOmega

    % =========================
    % 1) 输入检查 & 模式切换
    % =========================
    if numel(y) ~= 16, error('nondim_temp2 expects 16 inputs'); end
    state = y(1:15);

    if isempty(FixedOmega)
        W = y(16);          % 扫频
        current_Fw = Fw;
    else
        W = FixedOmega;     % 扫力
        current_Fw = y(16);
    end

    % =========================
    % 2) 解包参数
    % =========================
    be1 = sysP(1); be2 = sysP(2); mu_mass = sysP(3);
    al1 = sysP(4); ga1 = sysP(5); ze1 = sysP(6);
    
    lam_phys = sysP(7);   % 物理参数 lambda
    kap_e = sysP(8); kap_c = sysP(9); sigma = sysP(10); ga2 = sysP(11);

    % 计算对称耦合系数 theta = sqrt(lam)
    theta = sqrt(max(0, lam_phys)); 

    % =========================
    % 3) 状态向量拆分 & 基础算子
    % =========================
    x1 = state(1:5); x2 = state(6:10); q = state(11:15);
    W2 = W^2;
    
    x12 = x1 - x2;                    % 相对位移
    cubic12 = cubic_proj_013(x12);    % (x1-x2)^3
    cubic2  = cubic_proj_013(x2);     % x2^3

    % 二阶导算子 (3次谐波含9倍因子)
    M_op = [0; -W2; -W2; -9*W2; -9*W2];
    M_op_vec = @(u) M_op .* u;

    % =========================
    % 4) 构造耦合项向量
    % =========================
    % theta * q' (用于机械方程)
    theta_W = theta * W; 
    theta_W_3 = 3 * theta_W;
    
    % q' 的系数向量: [0, W*qs1, -W*qc1, 3W*qs3, -3W*qc3]
    force_em = [0; theta_W*q(3); -theta_W*q(2); theta_W_3*q(5); -theta_W_3*q(4)];

    % theta * (x1' - x2') (用于电路方程)
    x12_dot_coeff = [0; theta_W*x12(3); -theta_W*x12(2); theta_W_3*x12(5); -theta_W_3*x12(4)];

    % =========================
    % 5) 计算残差 R1, R2, R3
    % =========================
    
    % --- R1: 上层方程 ---
    % Residual = Inertia + Stiffness + Damping(EM) - ExtForce
    % 【修正】：这里必须是 + force_em (表现为正阻尼)
    R1 = M_op_vec(x1) + (be1+al1)*x12 + ga1*cubic12 + force_em;
    R1(2) = R1(2) - current_Fw; 

    % --- R2: 下层方程 ---
    % 阻尼项 2*mu*ze1*x2'
    c2_coef = 2*mu_mass*ze1*W;
    damp2 = [0; c2_coef*x2(3); -c2_coef*x2(2); 3*c2_coef*x2(5); -3*c2_coef*x2(4)];
    
    % 上层传递力 (注意 R2 中是 -UpperForce)
    Force_from_upper = (be1+al1)*x12 + ga1*cubic12;

    % 【修正】：这里必须是 - force_em (反作用力)
    R2 = mu_mass*M_op_vec(x2) + damp2 + be2*x2 + ga2*cubic2 ...
         - Force_from_upper - force_em;

    % --- R3: 电路方程 ---
    % Residual = Circuit_LHS - Source term
    % Source term = + theta*(x1'-x2')
    Circ_LHS = kap_e*M_op_vec(q) + kap_c*q;
    
    % 电阻项
    sig_W = sigma*W;
    Resist_term = [0; sig_W*q(3); -sig_W*q(2); 3*sig_W*q(5); -3*sig_W*q(4)];
    
    R3 = Circ_LHS + Resist_term - x12_dot_coeff;

    z = [R1; R2; R3];

    % =========================
    % 6) 解析雅可比 Jac
    % =========================
    if nargout > 1
        I5 = eye(5);
        Mat_Inertia = diag([0; -W2; -W2; -9*W2; -9*W2]);
        
        % 导数算子 D
        Mat_Deriv = zeros(5);
        Mat_Deriv(2,3) = W;    Mat_Deriv(3,2) = -W;
        Mat_Deriv(4,5) = 3*W;  Mat_Deriv(5,4) = -3*W;

        J_cubic_x12 = AFT_GetJac(x12);
        J_cubic_x2  = AFT_GetJac(x2);

        % --- R1 derivatives ---
        J11 = Mat_Inertia + (be1+al1)*I5 + ga1*J_cubic_x12;
        J12 = -(be1+al1)*I5 - ga1*J_cubic_x12;
        % 【修正】对应 R1 中的 +force_em，导数为 +theta*D
        J13 = theta * Mat_Deriv;

        % --- R2 derivatives ---
        J21 = J12; % Symmetric stiffness matrix
        J22 = mu_mass*Mat_Inertia + (2*mu_mass*ze1)*Mat_Deriv + be2*I5 + ga2*J_cubic_x2 ...
              + (be1+al1)*I5 + ga1*J_cubic_x12;
        % 【修正】对应 R2 中的 -force_em，导数为 -theta*D
        J23 = -theta * Mat_Deriv;

        % --- R3 derivatives ---
        % R3 ~ -theta * (x1' - x2')
        J31 = -theta * Mat_Deriv;
        J32 = +theta * Mat_Deriv;
        J33 = kap_e*Mat_Inertia + sigma*Mat_Deriv + kap_c*I5;

        Jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];
    end
end

% (AFT 辅助函数保持不变)
function cubic = cubic_proj_013(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    cubic  = T_inv * ((T_mat * u).^3);
end
function J_aft = AFT_GetJac(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    u_time = T_mat * u;
    df_du = 3 * u_time.^2;
    J_aft = T_inv * (df_du .* T_mat);
end
function [N, T_mat, T_inv] = get_AFT_matrices()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64; t = (0:pN-1)'*(2*pi/pN);
        c1=cos(t); s1=sin(t); c3=cos(3*t); s3=sin(3*t); dc=ones(pN,1);
        pT = [dc, c1, s1, c3, s3];
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv; pTinv(1,:) = (1/pN) * dc';
    end
    N = pN; T_mat = pT; T_inv = pTinv;
end