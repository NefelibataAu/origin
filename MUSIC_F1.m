function [P_MUSIC_dB,phi_e] = MUSIC_F1(k, K, d, X, P, L)
% MUSIC_F 
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：K:信号源的数量
% 输入：X(M*K):接收信号矩阵
% 输入：P:子阵数 L:子阵长度
% 输出：P_MUSIC_dB(1*len):谱估计函数  phi_e:信号源方向

[M, N] = size(X);
z = (0:d:(M-1)*d)'; 

%%% 阵列接收信号的协方差矩阵的特征分解
 R = X*X'/N;  
 R_f = zeros(size(L,L));
 for i = 1 : P
     X_temp = X(i:i+L-1, i:i+L-1);
     R_f = X_temp'*X_temp + R_f;
 end
 R_f = R_f./P;                          % 阵列接收信号的协方差矩阵
[EV, D] = eig(R_f);                     % 特征值分解
EVA = diag(D);                          % 提取特征值
[~, I] = sort(EVA, 'descend');          % 降序排序
Q = EV(:, I);                           % 特征向量构成的矩阵
Q_n = Q(:, K+1:L);   

%%% 计算MUSIC谱估计函数
phi_list = linspace(-pi/2, pi/2, 500)';
S1 = exp(-1j*k*z*sin(phi_list'));       % 不同方向对应的流型矢量构成矩阵
P_MUSIC = 1./sum(abs(Q_n'*S1).^2);      % MUSIC 谱估计公式

%%% 转换为dB
P_MUSIC = abs(P_MUSIC);
P_MUSIC_max = max(P_MUSIC);
P_MUSIC_dB = 10*log10(P_MUSIC/P_MUSIC_max);

%%% 提取信号源方向
[P_peaks, P_peaks_idx] = findpeaks(P_MUSIC_dB);     % 提取峰值
[~, I] = sort(P_peaks, 'descend');                  % 峰值降序排序
P_peaks_idx = P_peaks_idx(I);                       % 提取前M个
P_peaks_idx = P_peaks_idx(1:K);
phi_e = phi_list(P_peaks_idx)*180/pi;               % 估计方向
disp('MUSIC算法信号源估计方向为：');
disp(phi_e);


end

