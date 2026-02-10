function y=branch_aux2(x,sysP)

global tracking_file_name xc arc

% x 本身就是 continuation 的未知向量（应为 16 维）
% 关键：把完整 16 维向量传给 nondim_temp2
y1 = feval(tracking_file_name, x, sysP);

% 弧长约束也用完整向量
y = [y1; norm(x - xc) - arc];

end
