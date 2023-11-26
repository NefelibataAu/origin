function [P_MVDR_dB, phi_e] = MVDR_F(k, K, d, X)
% MVDR_F 
% 输入：k:信号波数 K:信号源的数量 d:接收阵列间隔 
% 输入：X(M*N):接收信号矩阵
% 输出：P_MVDR_dB(1*len):谱估计函数   phi_e:信号源方向

[M, N] = size(X);
z = (0:d:(M-1)*d)'; 


%% 算法估计
R = X*conj(X')/N; 
D = 500;
seita = linspace(-90, 90, D);
P_MVDR_dB = zeros(1, length(seita));
for i = 1 : length(seita)
    a_seita = exp( -1j*k*z*sind(seita(i)) );
    P_MVDR_dB(i) = 1/abs((conj(a_seita')*R^(-1)*a_seita));
end

%% 提取信号源方向
[P_peaks, P_peaks_idx] = findpeaks(P_MVDR_dB);      % 提取峰值
[~, I] = sort(P_peaks, 'descend');                  % 峰值降序排序
P_peaks_idx = P_peaks_idx(I);                       
P_peaks_idx = P_peaks_idx(1:K);
phi_e = seita(P_peaks_idx);                         % 估计方向
disp('MVDR算法信号源估计方向为：');
disp(phi_e');

end

