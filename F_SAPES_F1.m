function [P_SAPES_dB] = F_SAPES_F1(k, d, X, P, L, phi)

% F_SAPES得到波束图
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：phi:期望来波方向
% 输入：P:子阵数 L:子阵长度
% 输入：X(M*N):接收信号矩阵
% 输入：期望信号方向
% 输出：P_SAPES_dB(1*len):谱估计函数  

[M, N] = size(X);

%% 计算相关矩阵与平滑矩阵
 R = X*ctranspose(X)/N;   
 R_f = zeros(size(L,L));
 for i = 1 : P
     R_f = R(i:i+L-1,i:i+L-1) + R_f;
 end
 R_f = R_f./P; 
 
 %% 计算干扰相关矩阵
 z_p = ( 0:d:(P-1)*d )';
 a_p = ( exp(-1j*k*z_p*sin(phi')) )';
 T = zeros(M-P+1, M);
 for i = 1 : (M-P+1)
     T(i,:) = [zeros(1,i-1), a_p, zeros(1, M-P-i+1)];
 end
 G_f = T * R * ctranspose(T)/(P^(2));
 Q_f = R_f - G_f;
 
 %% 计算最优权向量
 z_L = (0:d:(L-1)*d)';
 a_L = exp(-1j*k*z_L*sin(phi'));
 omega = Q_f^(-1)*a_L/( ctranspose(a_L)*Q_f^(-1)*a_L );
 
 %% 得到波束图
D = 500;
seita = linspace(-90, 90, D);
F = abs( ctranspose(omega) *  exp(-1j*k*z_L*sind(seita)) );
P_SAPES_dB = 20*log10( F/(max (F) ) );

%%% 绘图
figure;
plot(seita, P_SAPES_dB);
xlabel('空间角度/(°)');
ylabel('归一化方向图/dB');
grid on;
hold on;

end

