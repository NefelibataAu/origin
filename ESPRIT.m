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
SNRs = [20, 20];      % 信噪比(dB)
N = 1000;             % 采样点数

%% 阵列接收信号仿真模拟
[X] = Signal_Generator(k, z, phi, SNRs, N);


[phi_e] = ESPRIT_F(k, K, d, X);



