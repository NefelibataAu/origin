function [P_SAPES_dB, phi_e] = F_SAPES_F(k, K, d, X, P, L)
% F_SAPES doa估计
% 输入：k:信号波数 K:信号源的数量 d:接收阵列间隔  
% 输入：X(M*N):接收信号矩阵
% 输入：P:子阵数 L:子阵长度
% 输出：P_SAPES:谱估计 phi_e:信号源方向 

[M, N] = size(X);
%% 计算平滑R_f
 R = X*X'/N;  
 R_f = zeros(size(L,L));
 for i = 1 : P
     X_temp = X(i:i+L-1, i:i+L-1);
     R_f = X_temp'*X_temp + R_f;
 end
 R_f = R_f./P;
 
 %% 谱估计计算
 D = 500;
 seita = linspace(-90, 90, D);
 P_SAPES_dB = zeros(1, D);
 z_p = ( 0:d:(P-1)*d )';
 for i = 1 : D
    a_p = ( exp(1j*k*z_p*sind(seita(i)')) )';
    T = zeros(M-P+1, M);
    for j = 1 : (M-P+1)
        T(j,:) = [zeros(1,j-1), a_p, zeros(1, M-P-j+1)];
    end
    G_f = T * R * transpose(conj(T))/(P^(2));
    Q_f = R_f - G_f;
    
    z_L = (0:d:(L-1)*d)';
    a_L = exp(-1j*k*z_L*sind(seita(i)'));
    omega = Q_f^(-1)*a_L/( ctranspose(a_L)*Q_f^(-1)*a_L );
    P_SAPES_dB(i) = abs( (conj(omega')*G_f*omega) );
    
 end
%% 提取信号源方向
[P_peaks, P_peaks_idx] = findpeaks(P_SAPES_dB);     % 提取峰值
[~, I] = sort(P_peaks, 'descend');                  % 峰值降序排序
P_peaks_idx = P_peaks_idx(I);                     
P_peaks_idx = P_peaks_idx(1:K);
phi_e = seita(P_peaks_idx);                         % 估计方向
disp('SAPES算法信号源估计方向为：');
disp(phi_e');
 

