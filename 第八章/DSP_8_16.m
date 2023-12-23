clear;
close all;
clc

%%  数据初始化
M = 10;                         % 阵元数
lamd = 2;                       % 波长
d = lamd / 2;                   % 单位间距 
Array = (0 : d : M - 1);        % 阵列结构
DOA = [-10 40];                 % 信号入射角
K = length(DOA);
SNR_dB = [10 20];               % 信噪比，单位:dB
N =100;                         % 快拍数
noise_power  = 1;               % 噪声功率
amp = sqrt(noise_power * 10 .^(SNR_dB / 10));      % 信号幅值     
A = zeros(M, K);              
step = 0.01;                    % 谱搜索的步长
DOA_grid = (-90 : step : 90);   % 角度网格
P_MUSIC = zeros(1, length(DOA_grid));
P_MVDR = zeros(1, length(DOA_grid));
DOA_RootMUSIC = zeros(1, K);
DOA_ESPRIT = zeros(1, K);

%%  构建接收信号矩阵
S = diag(amp) / sqrt(2) * ( randn(K, N) + 1i * randn(K, N) );  % 信号矩阵
V = sqrt(noise_power / 2) * ( randn(M, N) + 1i * randn(M, N));
for k = 1 : K                   % 计算方向矩阵
    A(:, k) = exp(-1i * (0 : M - 1) * 2 * pi * d * ...
    sin(DOA(k) * pi / 180) / lamd);
end
X = A * S + V;                   % 接收信号矩阵

%%  (1)MUSIC 算法估计DOA
R = X * X' / N;
[V, D] = eig(R);
[Y, I] = sort(diag(D), 'descend');
G = V(:, I(K + 1 : end));                       % 噪声子空间

for m = 1 : length(DOA_grid) 
    a = exp(-1i * (0 : M - 1) * 2 * pi * d * ...
        sin(DOA_grid(m) * pi / 180) / lamd).';
    P_MUSIC(m) = 1 / (a' * G * G' * a);
end

% 绘制MUSIC谱图
P_MUSIC =  abs(P_MUSIC) / max(abs(P_MUSIC));      % 归一化MUSIC谱
P_MUSIC_dB = 10 * log10(P_MUSIC);
figure;
plot(DOA_grid , P_MUSIC_dB);
set(gca, 'XTick', (-90 : 10 : 90));
xlim([-90 90]);
xlabel('DOA/degree');
ylabel('归一化MUSIC谱/dB');
title('MUSIC算法估计信号DOA');

%%  (2)Root-MUSIC 算法估计DOA
syms z;
a_zz = z .^(0 : M -1);
a_z = z .^(-(0 : M -1)).'; 
P_RootMUSIC = a_zz * G * G' * a_z;          % 构造多项式

z_root = roots(sym2poly(z .^(M - 1) * P_RootMUSIC));   % 求根
[t, Index] = sort(abs(abs(z_root) - 1));
for k = 1 : K
    DOA_RootMUSIC(k) = asin(angle(z_root(Index(2 * k - 1))) ...
        * lamd / (2 * pi * d)) * 180 / pi;
end

disp('Root-MUSIC算法, DOA:');
sort(DOA_RootMUSIC)

%%  (3)ESPRIT 算法估计DOA
S = V(:, I(1 : K));                    % 信号子空间
S1 = S(1 : M - 1, :);
S2 = S(2 : M, :);
fai = S1 \ S2;

[~, D_fai] = eig(fai);
D_fai = diag(D_fai);
for k = 1 : K
    DOA_ESPRIT(k) = asin(-angle(D_fai(k)) * lamd / (2 * pi * d)) * 180 / pi;
end

disp('ESPRIT算法, DOA:');
sort(DOA_ESPRIT)

%% (4)MVDR 算法估计DOA

for m = 1 : length(DOA_grid) 
    a = exp(-1i * (0 : M - 1) * 2 * pi * d * ...
        sin(DOA_grid(m) * pi / 180) / lamd).';
    P_MVDR(m) = 1 / (a' * inv(R) * a);
end

% 绘制MVDR谱图
P_MVDR =  abs(P_MVDR) / max(abs(P_MVDR));      % 归一化MUSIC谱
P_MVDR_dB = 10 * log10(P_MVDR);
figure;
plot(DOA_grid , P_MVDR_dB);
set(gca, 'XTick', (-90 : 10 : 90));
xlim([-90 90]);
xlabel('DOA/degree');
ylabel('归一化MVDR谱/dB');
title('MVDR算法估计信号DOA');

%% (5)F-SAPES 算法估计DOA
P = 6;                               
L = M + 1 - P;                       
Rf = zeros(L, L);
for i = 1 : P
    Rf = Rf + X(i : i + L -1, :) * X (i : i + L - 1, :)' / N;
end
Rf = Rf / P;                        
n1 = 0 : P - 1;
n2 = 0 : L - 1;
cc = [1 zeros(1, L - 1)];
for n3 = -90 : 0.5 : 90
    fy = exp(1i * pi * sin(n3 / 180 * pi));
    tt = [(fy .^(n1')).' zeros(1, M - P)];
    Tfy = toeplitz(cc, tt);            
    GfTheta = 1 ./ (P ^2) * Tfy * R * Tfy'; 
    Qf = Rf - GfTheta;                      
    aTheta = fy .^(-n2');
    Wof = ((inv(Qf)) * aTheta) ./ (aTheta' * inv(Qf) * aTheta); 
    sigma2sTheta(((n3 + 90) / .5 + 1)) = Wof' * GfTheta * Wof;
end

figure;
sigma2sTheta = abs(sigma2sTheta) / max(abs(sigma2sTheta));
sigma2sTheta_dB = 10 * log10(sigma2sTheta);
x_angle = linspace(-90 , 90, length(sigma2sTheta));
plot(x_angle, sigma2sTheta_dB);
xlim([-90 90]);
set(gca, 'XTick', -90 : 10 : 90);
xlabel('DOA/degree');
ylabel('归一化F-SAPES谱/dB');
title('F-SAPES算法估计信号DOA');

