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
P = 2;                  % 阵元分组
L = 15;                 % 每组阵元数
%%% 信号源与干扰源参数
sos = [-20]'*pi/180;        % 来波方向
Nos = length(sos);          % 信号源数目
soi = [60]'*pi/180;         % 干扰源来波
Noi = length(soi);          % 干扰源数目
phi = [sos; soi];
K = Nos + Noi;
%%% 仿真参数
SNRs = 20;               % 信噪比(dB)
SNRi = 40;
N = 500;                % 采样点数

%% 估计相关矩阵
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

%%% 计算相关矩阵与平滑矩阵
 R = X*X'/N;  
 R_f = zeros(size(L,L));
 for i = 1 : P
     X_temp = X(i:i+L-1, :);
     R_f = X_temp * X_temp' + R_f;
 end
 R_f = R_f./P;
 
 %% 计算干扰相关矩阵
 z_p = ( 0:d:(P-1)*d )';
 a_p = ( exp(1j*k*z_p*sin(sos')) )';
 T = zeros(M-P+1, M);
 for i = 1 : (M-P+1)
     T(i,:) = [zeros(1,i-1), a_p, zeros(1, M-P-i+1)];
 end
 G_f = T * R * transpose(conj(T))/(P^(2));
 Q_f = R_f - G_f;
 
 %% 计算最优权向量
 z_L = (0:d:(L-1)*d)';
 a_L = exp(-1j*k*z_L*sin(sos'));
 omega = Q_f^(-1)*a_L/( ctranspose(a_L)*Q_f^(-1)*a_L );
 
 %% 得到波束图
D = 500;
seita = linspace(-90, 90, D);
F =abs( ctranspose(omega) *  exp(-1j*k*z_L*sind(seita)) );
f = 10*log10( F/(max (F) ) );

%%% 绘图
figure;
plot(seita, f);
xlabel('空间角度/(°)');
ylabel('归一化方向图/dB');
grid on;
hold on;





