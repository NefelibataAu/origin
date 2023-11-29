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
L = 9;
%%% 信号源参数
phi1 = [-10]'*pi/180;      % 信号源来波
phi2 = [40]'* pi/180;
phi = [phi1 ; phi2];
K = length(phi);          % 信号源数目
%%% 仿真参数
SNR1 = 20;              % 信噪比(dB)
SNR2 = 40;
SNRs = [SNR1, SNR2];
N = 100;                % 采样点数

%% 阵列接收信号仿真模拟
[X] = Signal_Generator(k, z, phi, SNRs, N);

%% 采用RLS求最优权值









