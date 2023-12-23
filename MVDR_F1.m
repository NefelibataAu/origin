function [P_MVDR_dB] = MVDR_F1(k, d, phi, X)
% MVDR得到波束图：
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：phi:期望来波方向
% 输入：X(M*N):接收信号矩阵
% 输出：P_MVDR_dB(1*len):谱估计函数  

[M, N] = size(X);
z = (0:d:(M-1)*d)'; 

%% R与w估计
Aos = exp(-1j*k*z*sin(phi'));
R = X*ctranspose(X)/N;  
omega = R^(-1)*Aos/(ctranspose(Aos)*R^(-1)*Aos);

%% 计算归一化波束图
D = 500;
seita = linspace(-90, 90, D);
F = abs( ctranspose(omega) *  exp(-1j*k*z*sind(seita)) );
P_MVDR_dB = 20*log10( F/(max (F) ) );

%% 波束图
figure;
plot(seita, P_MVDR_dB);
xlabel('空间角度/(°)');
ylabel('归一化方向图/dB');
grid on;
hold on;

end

