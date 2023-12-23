clear;
close all;
clc

M = 10;
N = 100;
K = 2;
theta = [-10; 40] * pi /180;            % 期望角为-10，干扰为40
SNR = [10; 20];
Am = sqrt(10.^(SNR / 10));
S = Am * ones(1, N);
S(2, :) = S(2, :) .* exp(1i * 2 * pi * rand(1, N));
for a = 1 : M
    for b = 1 : K
        A(a, b) = exp(-1i * (a - 1) * pi * sin(theta(b)));
    end
end
V = zeros(M, N);
for m = 1 : M
    v = wgn(1, N, 0, 'complex');
    v = v -mean(v);
    v = v / std(v);
    V(m, :) = v;
end
X = A * S + V;



C = exp(-1i * (0 : M - 1)' * pi * sin(theta(1)));   % 约束矩阵
f = 1;                                              % 约束响应向量
wq = C * inv(C' * C) * f;
[m, n] = size(C);
[Q, R] = qr(C);                     % QR分解
Ca = Q(:, n + 1 : m);               % 约束矩阵的正交补矩阵
 
%% 利用RLS算法获得最优权向量
det = 0.004;
wa = zeros(M - 1, 1);                 % 初始化权向量
lamd = 0.98;
P = det * eye(M - 1);
for n = 1 : N
    u = Ca' * X(:, n);                 % 输入向量
    k = lamd^(-1) * P * u / (1 + lamd^(-1) *u' * P * u);
    d = wq' * X(:, n);                 % 期望响应
    epsilon = d - wa' * u;
    wa = wa + k * conj(epsilon);       % 权向量更新
    P = lamd^(-1) * P - lamd^(-1) * k * u' * P;
end
w = wq - Ca * wa;

%% 归一化波束图
y = [];
for n = -90 : 0.5 : 90
    a = exp(-1i * (0 : M - 1).' * pi * sin(n * pi / 180));
    y = [y, w' * a];
end
y = 20 * log10(abs(y) / max(abs(y)));
x = linspace(-90, 90, length(y));
figure;
plot(x, y);
xlim([-90 90]);
xlabel('DOA/degree');
ylabel('归一化方向图/dB');
title('RLS算法实现MVDR波束形成器');
    


