function [phi_e] = RootMUSIC_F(k, K, d, X)
% MUSIC_F 
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：K:信号源的数量
% 输入：X(M*K):接收信号矩阵
% 输出：phi_e:信号源方向

[M, N] = size(X);
z = (0:d:(M-1)*d)'; 

%%% 阵列接收信号的协方差矩阵的特征分解
R = X*X'/N;                             % 阵列接收信号的协方差矩阵
[EV, D] = eig(R);                       % 特征值分解
EVA = diag(D);                          % 提取特征值
[~, I] = sort(EVA, 'descend');          % 降序排序
Q = EV(:, I);                           % 特征向量构成的矩阵
Q_n = Q(:, K+1:M);   

%%% 求解高阶方程
syms z_root
a_z = z_root.^(-[1:M]');
equations = [conj(a_z')*Q_n*conj(Q_n')*a_z];
solutions = solve(equations, z_root);
solutions = double(solutions);

%%% 提取信号源方向
solution_change = abs(solutions);
solution_change = abs(solution_change-1);
[~, I] = sort(solution_change); 
phi = solutions(I);
phi = phi(1:K);
phi_e = -angle(phi)*180/(pi*k*d);

disp('Root_MUSIC算法信号源估计方向为：');
disp(phi_e');


end

