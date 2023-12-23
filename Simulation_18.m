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

%% LCMV
C = exp(-1j*k*z*sin(phi'));                % 约束
f = [1; 1];                                  % 约束结果
omega = C*inv(ctranspose(C)*C)*f;       
[m, n] = size(C);           
[Q, R] = qr(C);                           % QR分解        
Ca = Q(:, n + 1 : m);                     % 约束矩阵的正交补矩阵

%% 采用RLS求最优权值
det = 0.004;
lam = 0.98;
P = det * ones(M-2);
wa = zeros(M-2, 1);
for n = 1 : N
    u = Ca' * X(:, n);                      % 输入向量
    ki = lam^(-1) * P * u / (1 + lam^(-1) *u' * P * u);
    d = omega' * X(:, n);                   % 期望响应
    epsilon = d - wa' * u;
    wa = wa + ki * conj(epsilon);           % 权向量更新
    P = lam^(-1) * P - lam^(-1) * ki * u' * P;
end
omega = omega - Ca * wa;

D = 500;
seita = linspace(-pi/2, pi/2, D);
F = abs( ctranspose(omega) *  exp(-1j*k*z*sin(seita)) );
f = 10*log10( F/(max (F) ) );

%% 绘图
figure;
axis([-100 100 -60 0]);
plot(seita*180/pi, f);
xlabel('空间角度/(°)');
ylabel('归一化方向图/dB');
grid on;
hold on;








