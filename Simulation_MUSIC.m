clc; clear; close all;

%% 参数设置
%%% 信号参数
c = 3e8;
freq = 10e9;
lambda = c/freq;        % 波长
k = 2*pi/lambda;        % 波数

%%% 阵列参数
M = 8;                  % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(M-1)*d)';     % 阵元坐标分布
%%% 信号源参数
phi = [30, 0]'*pi/180;          
K = length(phi);        % 信号源数目
%%% 仿真参数
SNRs = [20, 20]';       % 信噪比(dB)
N = 1000;               % 快拍数

%% 阵列接收信号仿真模拟
[X] = Signal_Generator(k, z, phi, SNRs, N);

%% DOA估计
D = 500;
% FT_Space
P_FT_dB = FT_Space(k, z, X, D);
phi_FT_e = Search_phi(P_FT_dB, K);
disp('FT信号源估计方向为：');
disp(phi_FT_e);
draw_dB(P_FT_dB');
% MUSIC 算法
[P_MUSIC_dB] = MUSIC_F(k, K, d, X, D);
phi_MUSIC_e = Search_phi(P_MUSIC_dB, K);
disp('MUSIC信号源估计方向为：');
disp(phi_MUSIC_e);
draw_dB(P_MUSIC_dB');

