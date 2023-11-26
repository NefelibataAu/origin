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
phi = [10]'*pi/180;          
% 来波方向
N = length(phi);                % 信号源数目
%%% 仿真参数
SNR = 20;             % 信噪比(dB)
K = 1000;             % 采样点数

%% 阵列接收信号空间谱
A = exp(-1j*k*z*sin(phi'));         % 流型矩阵
S = randn(N, K);                    % 输入信号
X = A*S;                            % 阵列接收信号
P = fft(X);
%%% 转换为dB
P = abs(P);
P_max = max(P);
P_dB = 10*log10(P/P_max);   
%%% 绘图
figure;
fai = 180/M;
plot(P_dB, 'k', 'Linewidth', 2);
xlabel('\phi (deg)');
ylabel('Pseudo-spectrum (dB)');
axis([-90, 90, -40, 0]);





