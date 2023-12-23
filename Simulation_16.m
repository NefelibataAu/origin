
clc; clear; close all;

%% 参数设置
%%% 工作频率
c = 3e8;
freq = 10e9;
lambda = c/freq;    % 波长
k = 2*pi/lambda;    % 波数
%%% 阵列参数
M = 10;                 % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(M-1)*d)';     % 阵元坐标分布
P = 2;
L = M - P + 1;
%%% 信号源参数
phi1 = [-10]'*pi/180;      % 信号源来波
phi2 = [40]'* pi/180;
phi = [phi1 ; phi2];
K = length(phi);          % 信号源数目
%%% 仿真参数
SNR1 = 10;              % 信噪比(dB)
SNR2 = 20;
N = 100;                % 采样点数

%% 阵列接收信号仿真模拟
[X] = Signal_Generator(k, z, phi, [SNR1, SNR2], N);
D = 500;

%% 不同算法的DOA估计
[P_MUSIC_dB] = MUSIC_F(k, K, d, X, D);
phi_MUSIC_e = Search_phi(P_MUSIC_dB, K);
disp('MUSIC算法信号源估计方向为：');
disp(phi_MUSIC_e);

[phi_RootMUSIC_e] = RootMUSIC_F(k, K, d, X);

[phi_ESPRIT_e] = ESPRIT_F(k, K, d, X);

[P_MVDR_dB] = MVDR_doa(k, d, X, D);
phi_MVDR_e = Search_phi(P_MVDR_dB, K);
disp('MVDR算法信号源估计方向为：');
disp(phi_MVDR_e);

[P_SAPES_dB] = F_SAPES_doa(k, d, X, P, L, D);
phi_SAPES_e = Search_phi(P_SAPES_dB, K);
disp('SAPES算法信号源估计方向为：');
disp(phi_SAPES_e);
draw_dB([P_MUSIC_dB';P_MVDR_dB;P_SAPES_dB]);
legend( 'MUSIC','MVDR','F-SAPES')





