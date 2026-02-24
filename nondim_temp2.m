function [z, Jac] = nondim_temp2(y, sysP)
%% nondim_temp2：Two-stage（T1）无量纲 HB/AFT 残差与雅可比（物理一致版）
%
% 新增：当 lam_phys≈0（断开电磁耦合）时，加入层间阻尼 zeta12（作用于 x1'-x2'）
%      满足作用-反作用：R1 +Fd12，R2 -Fd12
%
% 输入:
%   y    : 16x1 = [x1(5); x2(5); q(5); Omega_or_Fw]
%   sysP : 11x1 参数向量
% 输出:
%   z    : 15x1 残差
%   Jac  : 15x15 雅可比

    global Fw FixedOmega

    %=========================
    % 1) 输入检查 & 模式切换
    %=========================
    if numel(y) ~= 16
        error('nondim_temp2 expects y to be 16x1: [x1(5); x2(5); q(5); Omega_or_Fw].');
    end
    if numel(sysP) ~= 11
        error('nondim_temp2 expects sysP to be 11x1.');
    end

    state = y(1:15);

    if isempty(FixedOmega)
        W = y(16);          % 扫频：y(16)=Omega
        current_Fw = Fw;    % 激励幅值固定
    else
        W = FixedOmega;     % 扫力：Omega固定
        current_Fw = y(16); % y(16)=Fw
    end

    %=========================
    % 2) 解包参数
    %=========================
    be1 = sysP(1);
    be2 = sysP(2);
    mu_mass = sysP(3);

    al1 = sysP(4);
    ga1 = sysP(5);
    ze1 = sysP(6);

    lam_phys = sysP(7);
    kap_e = sysP(8);
    kap_c = sysP(9);
    sigma = sysP(10);
    ga2 = sysP(11);

    % 对称耦合系数 theta = sqrt(lam)
    theta = sqrt(max(0, lam_phys));

    %=========================
    % 2.5) 新增：lam≈0 时层间阻尼
    %=========================
    lam_eps = 1e-12;
    zeta12 = 0.0;
    if abs(lam_phys) < lam_eps
        zeta12 = 0.05;  % 断耦合时的层间阻尼比（线性对比模型用）0.020412
    end

    %=========================
    % 3) 状态拆分 & 基础算子
    %=========================
    x1 = state(1:5);
    x2 = state(6:10);
    q  = state(11:15);

    x12 = x1 - x2;                  % 相对位移
    cubic12 = cubic_proj_013(x12);  % (x1-x2)^3 -> 投影到0/1/3
    cubic2  = cubic_proj_013(x2);   % x2^3      -> 投影到0/1/3

    W2 = W^2;

    % 二阶导数算子（系数域）：d2/dt2: cos/sin -> -W^2, 3次谐波 -> -9W^2
    Mat_Inertia = diag([0; -W2; -W2; -9*W2; -9*W2]);  % 5x5

    % 一阶导数算子（系数域）：对 [a0,a1,b1,a3,b3]
    Mat_Deriv = zeros(5);
    Mat_Deriv(2,3) = W;    Mat_Deriv(3,2) = -W;
    Mat_Deriv(4,5) = 3*W;  Mat_Deriv(5,4) = -3*W;

    I5 = eye(5);

    %=========================
    % 4) 耦合项（统一用算子表达）
    %=========================
    % q' = D q,  x12' = D x12
    q_dot   = Mat_Deriv * q;
    x12_dot = Mat_Deriv * x12;

    % 机械电磁力项：+theta*q' 作用于 m1；-theta*q' 作用于 m2
    force_em = theta * q_dot;

    % 电路耦合项：theta*(x1' - x2')（写到残差左边）
    x12_dot_coeff = theta * x12_dot;

    % lam≈0 时层间阻尼：Fd12 = 2*zeta12*(x1' - x2')
    damp12 = (2*zeta12) * x12_dot;  % 若zeta12=0，则自然为0向量

    % 下层对地阻尼：2*mu*ze1*x2'
    x2_dot = Mat_Deriv * x2;
    damp2  = (2*mu_mass*ze1) * x2_dot;

    %=========================
    % 5) 残差 R1, R2, R3
    %=========================

    % --- R1 上层 ---
    % x1'' + (be1+al1)*x12 + ga1*x12^3 + theta*q' + damp12 = Fw*cos
    R1 = (Mat_Inertia*x1) ...
         + (be1+al1)*x12 + ga1*cubic12 ...
         + force_em ...
         + damp12;

    % 激励只施加在 cos(Ωt) 的系数 a1（第2个分量）
    R1(2) = R1(2) - current_Fw;

    % --- R2 下层 ---
    % mu*x2'' + damp2 + be2*x2 + ga2*x2^3 - [(be1+al1)*x12 + ga1*x12^3] - theta*q' - damp12 = 0
    Force_from_upper = (be1+al1)*x12 + ga1*cubic12;

    R2 = mu_mass*(Mat_Inertia*x2) ...
         + damp2 ...
         + be2*x2 + ga2*cubic2 ...
         - Force_from_upper ...
         - force_em ...
         - damp12;

    % --- R3 电路 ---
    % 与论文式(24)一致： kap_e*q'' + sigma*q' + kap_c*q = -theta*(x1' - x2')
    % => 残差写成： kap_e*q'' + sigma*q' + kap_c*q + theta*(x1' - x2') = 0
    Circ_LHS   = kap_e*(Mat_Inertia*q) + kap_c*q;
    Resist_term = sigma * q_dot;

    R3 = Circ_LHS + Resist_term + x12_dot_coeff;

    z = [R1; R2; R3];

    %=========================
    % 6) 解析雅可比 Jac
    %=========================
    if nargout > 1
        J_cubic_x12 = AFT_GetJac(x12);
        J_cubic_x2  = AFT_GetJac(x2);

        % --- dR1/dx1, dR1/dx2, dR1/dq ---
        J11 = Mat_Inertia + (be1+al1)*I5 + ga1*J_cubic_x12;
        J12 = -(be1+al1)*I5 - ga1*J_cubic_x12;
        J13 = theta * Mat_Deriv;  % +theta*q'

        % 层间阻尼对 R1 的导数：+2*zeta12*D*(x1-x2)
        if zeta12 ~= 0
            J11 = J11 + (2*zeta12)*Mat_Deriv;
            J12 = J12 - (2*zeta12)*Mat_Deriv;
        end

        % --- dR2/dx1, dR2/dx2, dR2/dq ---
        J21 = -(be1+al1)*I5 - ga1*J_cubic_x12;

        J22 = mu_mass*Mat_Inertia ...
              + (2*mu_mass*ze1)*Mat_Deriv ...
              + be2*I5 + ga2*J_cubic_x2 ...
              + (be1+al1)*I5 + ga1*J_cubic_x12;

        J23 = -theta * Mat_Deriv;  % -theta*q'

        % 层间阻尼对 R2 的导数：-2*zeta12*D*(x1-x2)
        if zeta12 ~= 0
            J21 = J21 - (2*zeta12)*Mat_Deriv;
            J22 = J22 + (2*zeta12)*Mat_Deriv;
        end

        % --- dR3/dx1, dR3/dx2, dR3/dq ---
        % R3 = ... + theta*D*(x1-x2)
        J31 = +theta * Mat_Deriv;
        J32 = -theta * Mat_Deriv;
        J33 = kap_e*Mat_Inertia + sigma*Mat_Deriv + kap_c*I5;

        Jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];
    end
end

%=========================
% AFT 辅助函数
%=========================
function cubic = cubic_proj_013(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    cubic = T_inv * ((T_mat * u).^3);
end

function J_aft = AFT_GetJac(u)
    [~, T_mat, T_inv] = get_AFT_matrices();
    u_time = T_mat * u;
    df_du = 3 * u_time.^2;        % d(u^3)/du = 3u^2
    J_aft = T_inv * (df_du .* T_mat);
end

function [N, T_mat, T_inv] = get_AFT_matrices()
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64;
        t = (0:pN-1)'*(2*pi/pN);

        c1 = cos(t);  s1 = sin(t);
        c3 = cos(3*t); s3 = sin(3*t);
        dc = ones(pN,1);

        pT = [dc, c1, s1, c3, s3];

        % 与你原实现一致的投影逆
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;
        pTinv(1,:) = (1/pN) * dc';
    end
    N = pN; T_mat = pT; T_inv = pTinv;
end