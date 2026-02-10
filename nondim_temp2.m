function [z, Jac] = nondim_temp2(y, sysP)
%% nondim_temp2：T1（标准 Two-stage Raft）拓扑下的 HBM 残差 + 解析雅可比
%
% 【用途】
%   - 用谐波平衡法 (HBM, 仅保留 0 次 + 1 次 + 3 次谐波) 构造非线性代数方程 z=0
%   - 同时返回对“15 个状态系数”的解析 Jacobian (15x15)，用于 Newton / fsolve / 弧长延拓加速收敛
%
% 【输入】
%   y   : 16x1 向量
%         y(1:15)  = HBM 状态系数（x1, x2, q 各 5 个系数）
%         y(16)    = 延拓参数（扫频模式：Omega；扫力模式：Fw）
%
%         每个变量的 5 个 HBM 系数顺序统一为：
%         [DC, cos(1*Ωt), sin(1*Ωt), cos(3*Ωt), sin(3*Ωt)]
%         即：u(t) ≈ u0 + uc1*cos(Ωt) + us1*sin(Ωt) + uc3*cos(3Ωt) + us3*sin(3Ωt)
%
%   sysP: 11x1 参数向量（与你的论文/代码约定一致）
%         sysP = [be1, be2, mu, al1, ga1, ze1, lam, kap_e, kap_c, sigma, ga2]
%
% 【输出】
%   z   : 15x1 残差向量，按 [R1(1:5); R2(1:5); R3(1:5)] 叠加
%   Jac : 15x15 雅可比矩阵 dz/d(state)，state=[x1(5); x2(5); q(5)]
%
% 【物理模型（T1 标准 Two-stage）】
%   变量含义：
%     x1：上层质量（设备）绝对位移（无量纲）
%     x2：下层质量（筏/基础）绝对位移（无量纲）
%     q ：电荷变量（i = dq/dt），用于电路二阶形式（RLC/NIC）
%
%   线性+非线性耦合方程（无量纲形式，与你文档一致）：
%     (1) x1'' + (be1+al1)(x1-x2) + ga1(x1-x2)^3  - lam*q' = f*cos(Ωt)
%     (2) μ x2'' + 2μζ2 x2' + be2 x2 + ga2 x2^3  - [(be1+al1)(x1-x2)+ga1(x1-x2)^3] + lam*q' = 0
%     (3) kap_e q'' + sigma q' + kap_c q - (x1' - x2') = 0
%
%   其中：
%     - (be1+al1, ga1) 作用在相对位移 (x1-x2) —— 这就是 T1 “上层坐筏” 的关键
%     - (be2, ga2)     作用在 x2（下层对地）—— 下层对地线性/非线性支撑
%     - 电磁耦合项按审稿人要求：上层 -lam*q'，下层 +lam*q'，保证作用-反作用，且被动电阻 sigma>0 时为正阻尼
%
%   HBM 的做法：
%     把每个未知量展开成 0/1/3 次谐波系数，代入微分方程后，分别对每个基函数投影得到 15 个代数方程。
%
    global Fw FixedOmega

    % =========================
    % 1) 输入检查
    % =========================
    if numel(y) ~= 16, error('nondim_temp2 expects 16 inputs'); end
    if numel(sysP) ~= 11, error('sysP must be 11x1'); end

    state = y(1:15);

    % =========================
    % 2) 扫频/扫力模式切换
    % =========================
    % 约定：
    %   - FixedOmega 为空：执行扫频 FRF，y(16) 代表当前频率 Ω
    %   - FixedOmega 非空：执行定频扫力，y(16) 代表当前激励幅值 f(=Fw)
    if isempty(FixedOmega)
        W = y(16);          % 扫频：延拓参数为频率 Ω
        current_Fw = Fw;    % 激励幅值使用全局 Fw
    else
        W = FixedOmega;     % 扫力：频率固定为 FixedOmega
        current_Fw = y(16); % 延拓参数为激励幅值
    end

    % =========================
    % 3) 解包参数（顺序必须和 sysP 一致）
    % =========================
    be1 = sysP(1);  % 上下层之间的线性耦合刚度系数（通常归一化为 1）
    be2 = sysP(2);  % 下层对地线性刚度系数
    mu_mass = sysP(3); % 质量比 μ = m2/m1
    al1 = sysP(4);  % 上层（QZS）线性项系数（作用在 x1-x2 上）
    ga1 = sysP(5);  % 上层（QZS）立方非线性系数（作用在 x1-x2 上）
    ze1 = sysP(6);  % 下层对地阻尼比 ζ2（注意：实际方程是 2μζ2 x2'）
    lam = sysP(7);  % 机电耦合强度 λ（电磁力/反电动势归一化后系数）
    kap_e = sysP(8);% 电路“电感项”系数 κe
    kap_c = sysP(9);% 电路“倒电容项”系数 κc
    sigma = sysP(10);% 等效电阻项 σ（NIC 时可为负；被动电阻 σ>0）
    ga2 = sysP(11); % 下层对地立方非线性系数（作用在 x2 上）

    % =========================
    % 4) 状态向量拆分：x1, x2, q 各 5 个系数
    % =========================
    x1 = state(1:5);     % 上层位移 HBM 系数
    x2 = state(6:10);    % 下层位移 HBM 系数
    q  = state(11:15);   % 电荷 HBM 系数（电流 i = q'）

    W2 = W^2;

    % =========================================================================
    % PART 1：计算残差向量 z（15x1）
    % =========================================================================

    % =========================
    % 5) 一些常用中间量
    % =========================
    % 相对位移：x12 = x1 - x2（T1 的上层弹簧/非线性都作用在这里）
    x12 = x1 - x2;

    % 非线性项投影：把 (x12(t))^3 映射回 0/1/3 次谐波系数（AFT）
    cubic12 = cubic_proj_013(x12);

    % 下层对地非线性：x2^3 的 0/1/3 次谐波系数
    cubic2  = cubic_proj_013(x2);

    % =========================
    % 6) “二阶导数”在频域的对角算子（HBM 关键）
    % =========================
    % 若 u(t)=a cos(Ωt)+b sin(Ωt)，则 u''(t)=-(Ω^2)a cos(Ωt) -(Ω^2)b sin(Ωt)
    % 对 3 次谐波：频率是 3Ω，因此 (3Ω)^2 = 9Ω^2
    % 所以在系数空间：
    %   u'' -> [0, -Ω^2, -Ω^2, -9Ω^2, -9Ω^2] .* u
    % 【重要】这里必须包含 9 倍因子，否则残差与雅可比不一致，会导致牛顿/延拓不收敛或假解。
    M_op = [0; -W2; -W2; -9*W2; -9*W2];
    M_op_vec = @(u) M_op .* u;

    % =========================
    % 7) “一阶导数”在系数空间的显式形式（用于 q'、x'）
    % =========================
    % 对基函数 [DC, cos, sin, cos3, sin3]：
    %   d/dt(cos)= -Ω sin
    %   d/dt(sin)=  Ω cos
    % 因此：
    %   u' 的系数为：
    %     DC: 0
    %     cos1 系数来自 sin1： +Ω * us1
    %     sin1 系数来自 cos1： -Ω * uc1
    %     cos3 系数来自 sin3： +3Ω * us3
    %     sin3 系数来自 cos3： -3Ω * uc3
    %
    % 这里先组装 λ*q' 在各个谐波上的系数向量（后面 R1 用 -λq'，R2 用 +λq'）
    lam_W = lam*W; lam_W_3 = 3*lam_W;
    lam_dQ_1 = [0; lam_W*q(3); -lam_W*q(2); lam_W_3*q(5); -lam_W_3*q(4)];

    % =========================
    % 8) R1：上层方程残差（5x1）
    % =========================
    % 方程：x1'' + (be1+al1)(x1-x2) + ga1(x1-x2)^3 - lam*q' = f*cos(Ωt)
    %
    % 注意：激励力只作用在 cos(Ωt) 分量上，因此只在第 2 个系数位置扣除 current_Fw
    R1 = M_op_vec(x1) + (be1+al1)*x12 + ga1*cubic12 - lam_dQ_1;
    R1(2) = R1(2) - current_Fw;

    % =========================
    % 9) R2：下层方程残差（5x1）
    % =========================
    % 方程：μ x2'' + 2μζ2 x2' + be2 x2 + ga2 x2^3
    %        -[(be1+al1)(x1-x2)+ga1(x1-x2)^3] + lam*q' = 0
    %
    % 先构造阻尼项 2μζ2 x2' 在系数空间的贡献：
    c2_coeff = 2*mu_mass*ze1*W;
    c2_force = [0; c2_coeff*x2(3); -c2_coeff*x2(2); 3*c2_coeff*x2(5); -3*c2_coeff*x2(4)];

    % 上层弹簧对相对位移产生的力（在上层方程里是 +UpperForce，在下层方程里是 -UpperForce）
    UpperForce = (be1+al1)*x12 + ga1*cubic12;

    R2 = mu_mass*M_op_vec(x2) + c2_force + be2*x2 + ga2*cubic2 - UpperForce + lam_dQ_1;

    % =========================
    % 10) R3：电路方程残差（5x1）
    % =========================
    % 方程：kap_e*q'' + sigma*q' + kap_c*q - (x1' - x2') = 0
    %
    % 电路“弹性项”：kap_c*q
    Circ_LHS = kap_c*q;

    % 电路“电感项”：kap_e*q''（二阶导在频域使用 M_op，对 3 次谐波含 9 倍因子）
    Circ_LHS = Circ_LHS + kap_e*(M_op_vec(q));

    % 电路“电阻项”：sigma*q'
    sig_W = sigma*W;
    Circ_Damp = [0; sig_W*q(3); -sig_W*q(2); 3*sig_W*q(5); -3*sig_W*q(4)];

    % 机电耦合速度项：-(x1' - x2') = -d/dt(x1-x2)
    x12_dot = [0; W*x12(3); -W*x12(2); 3*W*x12(5); -3*W*x12(4)];

    R3 = Circ_LHS + Circ_Damp - x12_dot;

    % 总残差向量（15x1）
    z = [R1; R2; R3];

    % =========================================================================
    % PART 2：组装解析雅可比 Jac（15x15）
    % =========================================================================
    if nargout > 1
        % =========================
        % 11) 基本算子矩阵（5x5）
        % =========================
        I5 = eye(5);

        % 二阶导算子矩阵：对角阵 diag([0,-Ω^2,-Ω^2,-9Ω^2,-9Ω^2])
        Mat_Inertia = diag([0; -W2; -W2; -9*W2; -9*W2]);

        % 一阶导算子矩阵 D：u' = D*u
        Mat_Deriv = zeros(5);
        Mat_Deriv(2,3) = W;    Mat_Deriv(3,2) = -W;
        Mat_Deriv(4,5) = 3*W;  Mat_Deriv(5,4) = -3*W;

        % =========================
        % 12) 非线性项的雅可比（通过 AFT 计算）
        % =========================
        % 对 f(u)=u^3，有：df/du = 3u^2（在时域是对角阵）
        % 映射到系数域：J = T_inv * diag(3u(t)^2) * T_mat
        J_cubic_x12 = AFT_GetJac(x12); % (x1-x2)^3 对 (x1-x2) 的导数
        J_cubic_x2  = AFT_GetJac(x2);  % x2^3 对 x2 的导数

        % =========================
        % 13) 分块组装 Jac（每个块 5x5）
        % =========================

        % -------- R1 对 x1, x2, q 的偏导 --------
        % R1 = x1'' + (be1+al1)(x1-x2) + ga1*(x1-x2)^3 - lam*q' - f*cos
        J11 = Mat_Inertia + (be1+al1)*I5 + ga1*J_cubic_x12;     % dR1/dx1
        J12 = -(be1+al1)*I5 - ga1*J_cubic_x12;                  % dR1/dx2
        J13 = -lam * Mat_Deriv;                                 % dR1/dq

        % -------- R2 对 x1, x2, q 的偏导 --------
        % R2 = μ x2'' + 2μζ2 x2' + be2 x2 + ga2 x2^3 - UpperForce + lam*q'
        % UpperForce = (be1+al1)(x1-x2) + ga1(x1-x2)^3
        %
        % dR2/dx1 = - d(UpperForce)/dx1 = J12（与上面对称）
        J21 = J12;

        % dR2/dx2 = μ*Inertia + 2μζ2*Deriv + be2*I + ga2*J(x2^3)
        %           - d(UpperForce)/dx2
        % 由于 UpperForce 对 x2 的偏导是：-(be1+al1)I - ga1*J ；前面还有一个负号，因此变成 +(...)
        J22 = mu_mass*Mat_Inertia ...
              + (2*mu_mass*ze1) * Mat_Deriv ...
              + be2*I5 + ga2*J_cubic_x2 ...
              + (be1+al1)*I5 + ga1*J_cubic_x12;

        % dR2/dq = +lam*q' => +lam*D
        J23 = lam * Mat_Deriv;

        % -------- R3 对 x1, x2, q 的偏导 --------
        % R3 = kap_e*q'' + sigma*q' + kap_c*q - (x1' - x2')
        J31 = -Mat_Deriv;                                        % dR3/dx1
        J32 =  Mat_Deriv;                                        % dR3/dx2
        J33 = kap_e*Mat_Inertia + sigma*Mat_Deriv + kap_c*I5;     % dR3/dq

        % =========================
        % 14) 拼成最终 15x15 Jacobian
        % =========================
        Jac = [J11, J12, J13;
               J21, J22, J23;
               J31, J32, J33];
    end
end

% ==========================================================
% AFT（交替频时域法）辅助函数
% ==========================================================
function cubic = cubic_proj_013(u)
    % 【功能】给定系数 u（0/1/3 次谐波），计算 u(t)^3 的 0/1/3 次谐波系数
    %
    % 做法（AFT）：
    %   1) 用 T_mat 把系数域 -> 时域采样点：u_time = T_mat*u
    %   2) 在时域逐点做非线性：f_time = u_time.^3
    %   3) 用 T_inv 投影回系数域：cubic = T_inv*f_time
    [~, T_mat, T_inv] = get_AFT_matrices();
    u_time = T_mat * u;
    f_time = u_time.^3;
    cubic  = T_inv * f_time;
end

function J_aft = AFT_GetJac(u)
    % 【功能】返回 f(u)=u^3 在系数域的雅可比：d(coeff_out)/d(coeff_in)
    %
    % 理论：
    %   在时域：df/du = 3u(t)^2（对角阵）
    %   在系数域：J = T_inv * diag(3u(t)^2) * T_mat
    %
    % 实现优化：
    %   diag(df)*T_mat 等价于 (df .* T_mat)，逐行缩放即可，避免显式构造对角阵
    [~, T_mat, T_inv] = get_AFT_matrices();
    u_time = T_mat * u;            % N×1 时域采样
    df_du_time = 3 * u_time.^2;    % N×1

    J_aft = T_inv * (df_du_time .* T_mat); % (5×N)*(N×5) => 5×5
end

function [N, T_mat, T_inv] = get_AFT_matrices()
    % 【功能】构造并缓存 AFT 的“重构矩阵”和“投影矩阵”
    %
    % N       : 时域采样点数（越大越准，但计算更慢；一般 64 足够）
    % T_mat   : (N×5) 系数 -> 时域，u_time = T_mat * u_coeff
    % T_inv   : (5×N) 时域 -> 系数，u_coeff = T_inv * u_time
    %
    % 系数顺序：[DC, cos1, sin1, cos3, sin3]
    %
    persistent pN pT pTinv
    if isempty(pN)
        pN = 64;
        t = (0:pN-1)'*(2*pi/pN);

        c1 = cos(t);   s1 = sin(t);
        c3 = cos(3*t); s3 = sin(3*t);
        dc = ones(pN,1);

        % --- 重构矩阵：u(t)=DC + c1*C1 + s1*S1 + c3*C3 + s3*S3 ---
        pT = [dc, c1, s1, c3, s3];

        % --- 投影矩阵：把时域信号投影回 0/1/3 次谐波系数 ---
        % 对 DC：1/N * sum(u)
        % 对 cos/sin：2/N * sum(u*basis)
        Inv = [dc, 2*c1, 2*s1, 2*c3, 2*s3]';
        pTinv = (1/pN) * Inv;

        % 保险起见：显式确保 DC 行就是 1/N
        pTinv(1,:) = (1/pN) * dc';
    end

    N = pN; T_mat = pT; T_inv = pTinv;
end
