clc; clear; close all;

%% 参数设置
%%% 工作频率
c = 3e8;
freq = 10e9;
lambda = c/freq;    % 波长
k = 2*pi/lambda;    % 波数
%%% 阵列参数
M = 16;                 % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(M-1)*d)';     % 阵元坐标分布
P = 2;
L = M - P + 1 ;
%%% 信号源参数
phi1 = [-10]'*pi/180;       % 信号源来波
phi2 = [40]'* pi/180;
phi = [phi1 ; phi2];
K = length(phi);            % 信号源数目
%%% 仿真参数
SNR1 = 10;                  % 信噪比(dB)
SNR2 = 20;
N = 100;                    % 采样点数

%% 相干波信号接收模拟
A = exp(-1j*k*z*sin(phi1'));        % 流型矩阵
S = randn(length(phi1), N);         % 输入信号
X1 = A*S;                           % 阵列接收信号
X1 = awgn(X1, SNR1, 'measured');    % 加载高斯白噪声
A = exp(-1j*k*z*sin(phi2'));        % 流型矩阵
S = exp(1i*pi/6) * S;
X2 = A*S;                           % 阵列接收信号
X2 = awgn(X2, SNR2, 'measured');    % 加载高斯白噪声
X = X1 + X2;

%% DOA估计
D = 500;

[P_MUSIC_dB] = MUSIC_F(k, K, d, X, D);
[P_MUSICs_dB] = MUSIC_F1(k, K, d, X, P, L, D);
draw_dB([P_MUSICs_dB]');

[P_MVDR_dB] = MVDR_doa2(k, d, phi, X, P, L, D);

[P_SAPES_dB] = F_SAPES_doa(k, d, X, P, L, D);
draw_dB(P_SAPES_dB);










