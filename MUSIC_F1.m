function [P_MUSIC_dB] = MUSIC_F1(k, K, d, X, P, L, Point)
% MUSIC_F 
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：K:信号源的数量
% 输入：X(M*K):接收信号矩阵
% 输入：P:子阵数 L:子阵长度
% 输入：Point:谱估计点数
% 输出：P_MUSIC_dB(1*len):谱估计函数  

[~, N] = size(X);
z = (0:d:(L-1)*d)'; 

%%% 阵列接收信号的协方差矩阵的特征分解
 R = X*ctranspose(X)/N;      
 R_f = zeros(size(L,L));
 for i = 1 : P
     R_f = R(i:i+L-1,i:i+L-1) + R_f;
 end
 R_f = R_f./P;                          % 阵列接收信号的协方差矩阵
[EV, D] = eig(R_f);                     % 特征值分解
EVA = diag(D);                          % 提取特征值
[~, I] = sort(EVA, 'descend');          % 降序排序
Q = EV(:, I);                           % 特征向量构成的矩阵
Q_n = Q(:, K+1:L);   

%%% 计算MUSIC谱估计函数
seita = linspace(-pi/2, pi/2, Point)';
S1 = exp(-1j*k*z*sin(seita'));                                      % 不同方向对应的流型矢量构成矩阵
P_MUSIC = 1./(ctranspose(S1)*Q_n*ctranspose(Q_n)*S1);            % MUSIC 谱估计公式

%%% 转换为dB
P_MUSIC = abs(diag(P_MUSIC));
P_MUSIC_max = max(P_MUSIC);
P_MUSIC_dB = 10*log10(P_MUSIC/P_MUSIC_max);



end

