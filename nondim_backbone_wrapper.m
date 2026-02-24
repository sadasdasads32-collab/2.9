function R_bb = nondim_backbone_wrapper(y_bb, sysP)
    % 这是一个包装器，用于欺骗 branch_follow2，让它把振幅当作延拓参数
    % y_bb 的结构定义：
    % y_bb(1): x10
    % y_bb(2): Omega (原来这里是 a11，现在变成了未知数频率)
    % y_bb(3): b11
    % y_bb(4:15): 剩余的 HB 系数
    % y_bb(16): A (原来这里是 Omega，现在变成了延拓参数 a11)
    
    A     = y_bb(16);
    Omega = y_bb(2);
    
    % 1. 还原为原始的 16 维 y 向量，供给你的 nondim_temp2
    y_orig = zeros(16,1);
    y_orig(1)    = y_bb(1);
    y_orig(2)    = A;           % 将控制参数 A 赋给 a11
    y_orig(3)    = y_bb(3);
    y_orig(4:15) = y_bb(4:15);
    y_orig(16)   = Omega;       % 将未知的频率赋给原始的 Omega 位置
    
    % 2. 强制抹除所有耗散项和外激励 (骨架线必须是无阻尼自由振动)
    global Fw;
    Fw = 0;             % 外部激振力设为 0
    sysP_free = sysP;
    sysP_free(6)  = 0;  % ze1 (机械阻尼) 设为 0
    sysP_free(10) = 0;  % sigma (电路电阻) 设为 0
    
    % 3. 调用你原本的非线性残差方程
    R_orig = nondim_temp2(y_orig, sysP_free);
    
    % 4. 替换相位锚定方程
    % 当阻尼和外力为 0 时，x1 的正弦投影方程 (R_orig(3)) 自然满足。
    % 我们将其替换为 b11 = 0，从而锚定自由振动的相位，保证雅可比矩阵满秩。
    R_bb = R_orig;
    R_bb(3) = y_bb(3) - 0; 
end