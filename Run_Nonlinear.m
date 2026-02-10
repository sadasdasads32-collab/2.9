%% Run_Nonlinear_True_Check.m
% 审稿人修正版：激活双重非线性，确保看到弯曲
% 目的：克服“同向运动”导致的非线性隐身问题
clear all; clear functions; close all; clc; rehash;

clear functions; % 【关键】清除 persistent 变量，防止 AFT 维数冲突

% --- 1. 参数设置 (强非线性) ---
P.be1=1.0; P.be2=1.0; P.mu=0.5;
P.al1=0.0;   % 去掉线性 QZS，让非线性更纯粹
P.ga1=0.5;   % 【强非线性】上层
P.ze1=0.02;  % 【低阻尼】让共振峰更尖，非线性更容易被激发
P.lam=0.1;   % 弱机电耦合，避免电流阻尼掩盖非线性
P.kap_e=0.05; P.kap_c=1.0; P.sigma=0.5;
P.ga2=0.2;   % 【关键】开启下层对地非线性 (绝对位移项)

sysP = [P.be1, P.be2, P.mu, P.al1, P.ga1, P.ze1, P.lam, P.kap_e, P.kap_c, P.sigma, P.ga2];

global Fw FixedOmega
Fw = 0.05;       % 足够大的激振力
FixedOmega = []; 

% --- 2. 延拓求解器初始化 ---
w_start = 0.1;
y0 = zeros(15,1); 
opts = optimoptions('fsolve','Display','off','SpecifyObjectiveGradient',true);
[y_start(1:15), ~, flag] = fsolve(@(y) nondim_temp2([y; w_start], sysP), y0, opts);

if flag <= 0, error('Initial point failed!'); end
y_start(16) = w_start;

% --- 3. 弧长延拓循环 ---
ds = 0.05;          
ds_max = 0.1;       
MaxSteps = 800;     % 增加步数以覆盖高频

Result_List = [];   % [Omega, Amp_x1, Amp_x2, Amp_Rel]
Y_curr = y_start;   
Tangent = [zeros(15,1); 1]; 

fprintf('Running Nonlinear Check (Range 0.1 - 4.0)...\n');
fprintf('Looking for TWO resonance peaks with BENDING.\n');

for step = 1:MaxSteps
    % Predictor
    Y_pred = Y_curr + ds * Tangent;
    
    % Corrector
    [Y_new, converged] = Newton_Corrector(Y_pred, Tangent, sysP);
    
    if converged
        w_val = Y_new(16);
        
        % 计算幅值
        x1_vec = Y_new(1:5);
        x2_vec = Y_new(6:10);
        amp_x1 = norm(x1_vec(2:3)); 
        amp_x2 = norm(x2_vec(2:3));
        
        % 计算相对位移幅值 (x1-x2)
        x_rel_vec = x1_vec - x2_vec;
        amp_rel = norm(x_rel_vec(2:3));
        
        Result_List = [Result_List; w_val, amp_x1, amp_x2, amp_rel];
        
        % Update
        diff_vec = Y_new - Y_curr;
        Tangent = diff_vec / norm(diff_vec);
        Y_curr = Y_new;
        
        % Step control
        if step > 5, ds = min(ds * 1.1, ds_max); end
        
        if w_val > 4.0, break; end % 扫到 4.0
    else
        ds = ds * 0.5;
        if ds < 1e-4, break; end
    end
end

% --- 4. 绘图 (带解释) ---
figure('Position',[100,100,1200,500]);

% 子图1: X1 幅频响应
subplot(1,3,1);
plot(Result_List(:,1), Result_List(:,2), 'r.-', 'LineWidth', 1.5);
xlabel('\Omega'); ylabel('|X1|'); title('X1: Frequency Hardening');
grid on; axis tight;
% 标注两个峰
text(0.8, max(Result_List(:,2))*0.2, 'Mode 1 (In-Phase)', 'Color','b');
text(2.0, max(Result_List(:,2))*0.2, 'Mode 2 (Anti-Phase)', 'Color','b');

% 子图2: 相对位移 (揭示非线性源)
subplot(1,3,2);
plot(Result_List(:,1), Result_List(:,4), 'k.-', 'LineWidth', 1.5);
xlabel('\Omega'); ylabel('|X1 - X2|'); title('Relative Disp (Activates \gamma_1)');
grid on; axis tight;

% 子图3: 相位差 (辅助验证)
subplot(1,3,3);
plot(Result_List(:,1), Result_List(:,2), 'r--', 'LineWidth',1); hold on;
plot(Result_List(:,1), Result_List(:,3), 'b--', 'LineWidth',1);
legend('|X1|', '|X2|');
title('Comparison X1 vs X2'); grid on; axis tight;

% =========================================================================
% Newton Corrector (Same as before)
% =========================================================================
function [Y_sol, success] = Newton_Corrector(Y_pred, Tangent, sysP)
    tol = 1e-6; 
    max_iter = 12;
    success = false;

    % ====== 【强制：所有输入都变成 16×1】======
    Y_pred0 = Y_pred(:);
    T       = Tangent(:);

    if numel(Y_pred0) < 16
        error('Y_pred has less than 16 elements. It is corrupted.');
    end
    if numel(T) < 16
        error('Tangent has less than 16 elements. It is corrupted.');
    end

    Y_pred0 = Y_pred0(1:16);   % 裁剪为16
    T       = T(1:16);

    Y = Y_pred0;

    % 切向量归一化，防止数值炸
    nT = norm(T);
    if nT < 1e-14
        T = [zeros(15,1); 1];
        nT = norm(T);
    end
    T = T / nT;

    for iter = 1:max_iter
        % ====== 再次确保 Y 是 16×1 ======
        Y = Y(:);
        if numel(Y) < 16
            error('Y has less than 16 elements during Newton iterations.');
        end
        Y = Y(1:16);

        % ====== 现在拼接永远不会 vertcat 出错 ======
        x = Y(1:15);
        w = Y(16);

        [F, Jx] = nondim_temp2([x; w], sysP);   % F:15×1, Jx:15×15

        % 弧长约束：reference 固定在 predictor 点
        g = dot(T, (Y - Y_pred0));              % 标量
        Phi = [F; g];                            % 16×1

        if norm(Phi) < tol
            success = true;
            Y_sol = Y;
            return;
        end

        % dF/dw：有限差分
        delta_w = 1e-7 * max(1, abs(w));
        Fp = nondim_temp2([x; w + delta_w], sysP); % 单输出
        Jw = (Fp - F) / delta_w;                   % 15×1

        Jac_Aug = [Jx, Jw; T'];                    % 16×16

        dY = Jac_Aug \ Phi;
        Y  = Y - dY;
    end

    Y_sol = Y;
end
