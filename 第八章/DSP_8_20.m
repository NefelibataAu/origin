clear;
close all;
clc

%%  数据初始化
M = 8;                         % 阵元数
lamd = 2;                       % 波长
d = lamd / 2;                   % 单位间距 
Array = (0 : d : M - 1);        % 阵列结构
DOA = 10;                 % 信号入射角
K = length(DOA);
SNR_dB = 20;               % 信噪比，单位:dB
N =1000;                         % 快拍数
noise_power  = 1;               % 噪声功率
amp = sqrt(noise_power * 10 .^(SNR_dB / 10));      % 信号幅值     
A = zeros(M, K);              
step = 0.01;                    % 谱搜索的步长
DOA_grid = (-90 : step : 90);   % 角度网格
P_MUSIC = zeros(1, length(DOA_grid));
W = zeros(M, M);

%%  构建接收信号矩阵
S = diag(amp) / sqrt(2) * ( randn(K, N) + 1i * randn(K, N) );  % 信号矩阵
V = sqrt(noise_power / 2) * ( randn(M, N) + 1i * randn(M, N));
for k = 1 : K                   % 计算方向矩阵
    A(:, k) = exp(-1i * (0 : M - 1) * 2 * pi * d * ...
    sin(DOA(k) * pi / 180) / lamd);
end
X = A * S + V;                   % 接收信号矩阵

%% 归一化波束形成矩阵
B = M;
m = 0;
aa = exp(- 1i * pi);
for k = 0 : M - 1
    W(:, k + 1) = aa .^((0 : M - 1) * k * (2 / M));
end
T = 1 / sqrt(M) * W(:, m + 1 : m + B);

%%  在波束空间应用MUSIC算法实现DOA估计
y = T' * X;
R = y * y' / N;
[V, D] = eig(R);
[Y, I] = sort(diag(D));
G = V(:, I(B - K : -1 : 1));

for m = 1 : length(DOA_grid) 
    a1 = exp(-1i * (0 : M - 1) * 2 * pi * d * ...
        sin(DOA_grid(m) * pi / 180) / lamd).';
    a2 = T' * a1;
    P_MUSIC(m) = 1 / (a2' * G * G' * a2);
end

%% 绘制MUSIC谱图
P_MUSIC =  abs(P_MUSIC) / max(abs(P_MUSIC));      % 归一化MUSIC谱
P_MUSIC_dB = 10 * log10(P_MUSIC);
figure;
plot(DOA_grid , P_MUSIC_dB);
set(gca, 'XTick', (-90 : 10 : 90));
xlim([-90 90]);
xlabel('DOA/degree');
ylabel('归一化MUSIC谱/dB');
title('在波束空间应用MUSIC算法实现DOA估计');

