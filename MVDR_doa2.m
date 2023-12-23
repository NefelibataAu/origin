function [P_MVDR_dB] = MVDR_doa2(k, d, phi, X, P, L, Point)
% MVDR得到波束图：
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：K:信号源的数量
% 输入：X(M*N):接收信号矩阵
% 输入：P:子阵数 L:子阵长度
% 输入：Point:谱估计点数
% 输出：P_MVDR_dB(1*len):谱估计函数  

[M, N] = size(X);
z = (0:d:(L-1)*d)'; 

%% 功率计算
R = X*ctranspose(X)/N;      
R_f = zeros(size(L,L));
for i = 1 : P
    R_f = R(i:i+L-1,i:i+L-1) + R_f;
end
R_f = R_f./P;      
seita = linspace(-90, 90, Point);
P_MVDR_dB = zeros(1, length(seita));
for i = 1 : length(seita)
    a_seita = exp( -1j*k*z*sind(seita(i)) );
    P_MVDR_dB(i) = 1/abs((ctranspose(a_seita)*R_f^(-1)*a_seita));
end


%% 波束图
figure;
plot(seita, P_MVDR_dB);
xlabel('空间角度/(°)');
ylabel('归一化方向图/dB');
grid on;
hold on;

end

