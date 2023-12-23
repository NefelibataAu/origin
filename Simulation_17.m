clc; clear; close all;

%% 参数设置
%%% 工作频率
c = 3e8;
freq = 10e9;
lambda = c/freq;    % 波长
k = 2*pi/lambda;    % 波数
%%% 阵列参数
M = 8;                 % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(M-1)*d)';     % 阵元坐标分布
P = 2;
L = 7;
%%% 信号源与干扰源参数
sos = [0]'*pi/180;         % 来波方向
Nos = length(sos);         % 信号源数目
soi = [-60, 50]'*pi/180;   % 干扰源来波
Noi = length(soi);         % 干扰源数目
phi = [sos; soi];
K = Nos + Noi;
%%% 仿真参数
SNRs = 0;             % 信噪比(dB)
SNRi1 = 40;
SNRi2 = 20;
N = 100;              % 采样点数

%% 接收模拟
[X] = Signal_Generator(k, z, phi, [SNRs, SNRi1, SNRi2], N);

P_MVDR_dB = MVDR_F1(k, d, sos, X);
P_SAPES_dB = F_SAPES_F1(k, d, X, P, L, sos);


