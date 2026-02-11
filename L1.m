function x = L1(x_interest, sysP, varargin)
%% Level-1 Continuation: fixed Omega, sweep Force (S-curve)
% Input:
%   x_interest : can be 15x1 (coeffs only) OR 16x1 ([coeffs; Omega])
%                can also be row vector, will be reshaped
%   sysP       : parameter vector passed to nondim_temp2
% Optional:
%   varargin{1}: if abs(val) > F_switch => treat as init Force
%               else treat as step size
%   varargin{2}: branch_len (max continuation steps)
%
% Output:
%   x : 16xM solution branch, last row is Force (Fw)

    global Fw FixedOmega

    %% -----------------------------
    % 0) Robust defaults & thresholds
    % -----------------------------
    step_default = 1e-3;
    step_min     = 1e-6;
    F_switch     = 0.1;      % abs(val) > F_switch => init force
    branch_len_default = 2000;

    current_F = Fw;
    step      = step_default;

    %% -----------------------------
    % 1) Parse optional arguments
    % -----------------------------
    if nargin >= 3 && ~isempty(varargin)
        input_val = varargin{1};
        if ~isscalar(input_val) || ~isfinite(input_val)
            error('L1: varargin{1} must be a finite scalar.');
        end

        if abs(input_val) > F_switch
            current_F = input_val;   % treat as initial force
        else
            step = input_val;        % treat as step
        end
    end

    if step == 0
        error('L1: step cannot be 0.');
    end
    if abs(step) < step_min
        step = sign(step) * step_min;    % allow negative (reverse sweep) if needed
    end

    % --- NEW: branch_len optional ---
    branch_len = branch_len_default;
    if nargin >= 4 && ~isempty(varargin{2})
        bl = varargin{2};
        if ~isscalar(bl) || ~isfinite(bl) || bl < 50
            error('L1: branch_len must be a finite scalar >= 50.');
        end
        branch_len = round(bl);
    end

    %% -----------------------------
    % 2) Normalize input shape
    % -----------------------------
    x_interest = x_interest(:);   % force column

    if numel(x_interest) < 15
        error('L1: x_interest must have at least 15 elements (coeffs). Got %d.', numel(x_interest));
    end

    % Always take coeffs as 15x1 column (DO NOT transpose)
    x_coeffs = x_interest(1:15);

    %% -----------------------------
    % 3) Determine fixed Omega safely
    % -----------------------------
    fixed_omega = [];

    if numel(x_interest) >= 16
        fixed_omega = x_interest(16);   % standard: [coeffs; Omega]
    elseif numel(x_interest) >= 31
        fixed_omega = x_interest(31);   % legacy hook
    else
        if isempty(FixedOmega)
            error(['L1: Omega is missing. Provide x_interest=[coeffs;Omega] (16x1) ' ...
                   'or set global FixedOmega before calling L1.']);
        end
        fixed_omega = FixedOmega;
    end

    % Lock frequency for nondim_temp2 "fixed omega" mode
    FixedOmega = fixed_omega;

    %% -----------------------------
    % 4) Build start vector for Newton
    % -----------------------------
    % y_start is 16x1: [15 coeffs; Force]
    y_start = [x_coeffs; current_F];

    %% -----------------------------
    % 5) Newton correction (optional but helpful)
    % -----------------------------
    try
        x0 = newton('nondim_temp2', y_start, sysP);
    catch ME
        fprintf('L1: Newton correction failed (%s). Start from raw guess.\n', ME.message);
        x0 = y_start;
    end

    mu0 = x0(end);
    mu1 = mu0 + step;

    % Slightly perturb the second initial point to avoid degenerate tangent
    x1 = x0(1:15);
    eps_dir = 1e-6;
    x1 = x1 + eps_dir * randn(size(x1));  %#ok<RAND>

    fprintf('L1: start | Omega=%.6f | F0=%.6f -> F1=%.6f | step=%.2e | branch_len=%d\n', ...
            fixed_omega, mu0, mu1, step, branch_len);

    %% -----------------------------
    % 6) Continuation
    % -----------------------------
    [x, ~] = branch_follow2('nondim_temp2', branch_len, mu0, mu1, x0(1:15), x1, sysP);

    fprintf('L1: finished | steps=%d | F_end=%.6f\n', size(x,2), x(end,end));

    %% -----------------------------
    % 7) Clear global omega lock
    % -----------------------------
    FixedOmega = [];
end
