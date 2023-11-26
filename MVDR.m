clc; clear; close all;

%% 参数设置
%%% 工作频率
c = 3e8;
freq = 10e9;
lambda = c/freq;    % 波长
k = 2*pi/lambda;    % 波数
%%% 阵列参数
M = 16;                 % 阵元数量
d = 0.5*lambda;         % 阵元间隔 
z = (0:d:(M-1)*d)';     % 阵元坐标分布
%%% 信号源与干扰源参数
sos = [20]'*pi/180;        % 来波方向
Nos = length(sos);         % 信号源数目
soi = [-40, 60]'*pi/180;   % 干扰源来波
Noi = length(soi);         % 干扰源数目
phi = [sos; soi];
K = Nos + Noi;
%%% 仿真参数
SNRs = 0;             % 信噪比(dB)
SNRi = 40;
N = 500;             % 采样点数

%% 估计相关矩阵与权向量
%%% 接收情况
Aos = exp(-1j*k*z*sin(sos'));       % 流型矩阵
S = randn(Nos, N);                  % 输入信号
X = Aos*S;                          % 阵列接收信号
X1 = awgn(X, SNRs, 'measured');     % 加载高斯白噪声

Aoi = exp(-1j*k*z*sin(soi'));       % 流型矩阵
S = randn(Noi, N);                  % 输入信号
X = Aoi*S;                          % 阵列接收信号
X2 = awgn(X, SNRi, 'measured');     % 加载高斯白噪声

X = X1 + X2;

%%% R与w估计
R = X*X'/N;  
omega = R^(-1)*Aos/(ctranspose(Aos)*R^(-1)*Aos);

%% 计算归一化波束图
D = 500;
seita = linspace(-90, 90, D);
F =abs( ctranspose(omega) *  exp(-1j*k*z*sind(seita)) );
f = 20*log10( F/(max (F) ) );

%%% 绘图
figure;
plot(seita, f);
xlabel('空间角度/(°)');
ylabel('归一化方向图/dB');
grid on;
hold on;




