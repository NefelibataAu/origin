clc; clear; close all;

%% 参数设置
%%% 工作频率
c = 3e8;
freq = 10e9;
lambda = c/freq;    % 波长
k = 2*pi/lambda;    % 波数
%%% 阵列参数
M = 8;                  % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(M-1)*d)';     % 阵元坐标分布
%%% 信号源参数
phi = [0, 10]'*pi/180;     % 来波方向
K = length(phi);           % 信号源数目
%%% 仿真参数
SNR = 20;             % 信噪比(dB)
N = 1000;             % 采样点数

%% 蒙塔卡罗
seita = zeros(1,K);
num = 100;
for i = 1:num
   %%% 阵列接受信号仿真模拟
   A = exp(-1j*k*z*sin(phi'));     
   S = randn(K, N);
   X = A * S;
   X1 = awgn(X, SNR, 'measured');
   R = X1*X1'/N;                           % 阵列接收信号的协方差矩阵
   [EV, D] = eig(R);                       % 特征值分解
   EVA = diag(D);                          % 提取特征值
   [EVA, I] = sort(EVA, 'descend');        % 降序排序
   Q = EV(:, I);                           % 特征向量构成的矩阵
   sign = Q(:, 1:K);                       % 信号子空间
 
   %%% 构造并分块
   S1 = sign(1:(M-1), :);                  % Sign的前(M-1)行
   S2 = sign(2:M, :);                      % Sign的后(M-1)行
   S12 = [S1 S2];                          % 构造S12
   [U, V] = eig(conj(S12.')*S12);          % 特征值分解
   U11 = U(1:K, 1:K);
   U12 = U(1:K, K+1:2*K);
   U21 = U(K+1:2*K, 1:K);
   U22 = U(K+1:2*K, K+1:2*K);
    
  %%% 对Fai矩阵特征值分解得解
   Fai = U12 * (U22^(-1));
   [EV, D] = eig(Fai);
   seita_tmp = [angle(D(1)), angle(D(4))]*180/(pi*k*d);
   seita_tmp = sort(seita_tmp);
   seita = (-seita_tmp + seita);
end
seita = seita/num



