function [P_MVDR_dB] = MVDR_doa(k, d, X, Point)
% MVDR_F 
% 输入：k:信号波数 K:信号源的数量 d:接收阵列间隔 
% 输入：X(M*N):接收信号矩阵
% 输入：Point:谱估计点数
% 输出：P_MVDR_dB(1*len):谱估计函数   phi_e:信号源方向

[M, N] = size(X);
z = (0:d:(M-1)*d)'; 


%% 算法估计
R = X*ctranspose(X)/N; 
seita = linspace(-pi/2, pi/2, Point);
P_MVDR_dB = zeros(1, length(seita));
for i = 1 : length(seita)
    a_seita = exp( -1j*k*z*sin(seita(i)) );
    P_MVDR_dB(i) = 1/abs((ctranspose(a_seita)*R^(-1)*a_seita));
end

end

