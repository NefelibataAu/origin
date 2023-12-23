function [P_SAPES_dB] = F_SAPES_doa(k, d, X, P, L, Point)
% F_SAPES doa估计
% 输入：k:信号波数 K:信号源的数量 d:接收阵列间隔  
% 输入：X(M*N):接收信号矩阵
% 输入：P:子阵数 L:子阵长度
% 输入：Point:谱估计点数
% 输出：P_SAPES:谱估计 

[M, N] = size(X);
%% 计算平滑R_f
 R = X*ctranspose(X)/N;      
 R_f = zeros(size(L,L));
 for i = 1 : P
     R_f = R(i:i+L-1,i:i+L-1) + R_f;
 end
 R_f = R_f./P;                          % 阵列接收信号的协方差矩阵
     
 %% 谱估计计算
 seita = linspace(-90, 90, Point);
 P_SAPES_dB = zeros(1, Point);
 z_p = ( 0:d:(P-1)*d )';
 for i = 1 : Point
    a_p = ( exp(-1j*k*z_p*sind(seita(i)')) )';
    T = zeros(M-P+1, M);
    for j = 1 : (M-P+1)
        T(j,:) = [zeros(1,j-1), a_p, zeros(1, M-P-j+1)];
    end
    G_f = T * R * ctranspose(T)/(P^(2));
    Q_f = R_f - G_f;
    
    z_L = (0:d:(L-1)*d)';
    a_L = exp(-1j*k*z_L*sind(seita(i)'));
    omega = Q_f^(-1)*a_L/( ctranspose(a_L)*Q_f^(-1)*a_L );
    P_SAPES_dB(i) = abs( (ctranspose(omega)*G_f*omega) );
 end
 
end

 

