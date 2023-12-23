clear; 
close all;
clc;

%% 数据初始化
M = 8;                                    % 阵元数
N = 100;                                  % 样本数
NN = 1000;
theta_grid = linspace(-90, 90, NN);       % 遍历的角度
lamd = 2;                                 % 波长
d = lamd / 2;                             % 单位间距
Array = (0 : M -1) * d;                   % 阵列结构
F = zeros(1, NN);                         % 初始化方向图
 
theta = 0;
theta_dst1 = -60;
theta_dst2 = 50;
SNR = 0;
SNR_dst1 = 40;
SNR_dst2 = 20;
f = 0.15;
f_dst1 = 0.1;
f_dst2 = 0.2;

%%  期望信号与干扰信号
sigma_square = 1;
noise = sqrt(sigma_square / 2) * (randn(M, N) + 1i *randn(M, N));
amp = sqrt(sigma_square *10^(SNR / 10));
amp_dst1 = sqrt(sigma_square *10^(SNR_dst1 / 10));
amp_dst2 = sqrt(sigma_square *10^(SNR_dst2 / 10));

signal = amp * exp(1i * 2 * pi * f * (0 : N - 1) + 1i * 2 * pi * rand);   % 期望信号
signal_dst1 = amp_dst1 * exp(1i * 2 * pi * f_dst1 * (0 : N - 1) + 1i * 2 * pi * rand); % 干扰信号1
signal_dst2 = amp_dst2 * exp(1i * 2 * pi * f_dst2 * (0 : N - 1) + 1i * 2 * pi * rand); % 干扰信号2

%%  阵列接收信号
a = exp(-1i * 2 * pi * (0 : M - 1) *d * sin(theta * pi / 180) /lamd).';
a_dst1 = exp(-1i * 2 * pi * (0 : M - 1) *d * sin(theta_dst1 * pi / 180) /lamd).';
a_dst2 = exp(-1i * 2 * pi * (0 : M - 1) *d * sin(theta_dst2 * pi / 180) /lamd).';

x = a * signal + a_dst1 * signal_dst1 + a_dst2 * signal_dst2 + noise;

%% MVDR 算法估计最优权向量
R = x * x' /N;
w = inv(R) * a / (a' * inv(R) * a);

%% 计算方向图及绘图
for k = 1 : NN
    a_grid = exp(-1i * 2 * pi * (0 : M - 1) * d * sin(theta_grid(k) * pi / 180) / lamd).';
    F(k) = abs(w' * a_grid);
end

F = F / max(F);                       % 归一化方向图
F_dB = 20 * log10(F);                 % 取对数

figure;
plot(theta_grid, F_dB);
xlim([-90 90]);
xlabel('DOA/degree');
ylabel('归一化方向图/dB');
title('MVDR算法');

