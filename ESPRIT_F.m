function [phi_e] = ESPRIT_F(k, K, d, X)
% ESPRIT_F 
% 输入：k:信号波数 d:接收阵列间隔  
% 输入：K:信号源的数量
% 输入：X(M*N):接收信号矩阵
% 输出：phi_e:信号源方向

[M, N] = size(X);
z = (0:d:(M-1)*d)'; 

%% 蒙特卡洛
seita = zeros(1,K);
num = 100;
for i = 1:num
   %%% 阵列接受信号仿真模拟
   R = X*X'/N;                             % 阵列接收信号的协方差矩阵
   [EV, D] = eig(R);                       % 特征值分解
   EVA = diag(D);                          % 提取特征值
   [~, I] = sort(EVA, 'descend');          % 降序排序
   Q = EV(:, I);                           % 特征向量构成的矩阵
   sign = Q(:, 1:K);                       % 信号子空间
 
   %%% 构造并分块
   S1 = sign(1:(M-1), :);                  % Sign的前(M-1)行
   S2 = sign(2:M, :);                      % Sign的后(M-1)行
   S12 = [S1 S2];                          % 构造S12
   [U, ~] = eig(conj(S12.')*S12);          % 特征值分解
   U12 = U(1:K, K+1:2*K);
   U22 = U(K+1:2*K, K+1:2*K);
    
  %%% 对Fai矩阵特征值分解得解
   Fai = U12 * (U22^(-1));
   [~, D] = eig(Fai);
   seita_tmp = [angle(D(1)), angle(D(4))]*180/(pi*k*d);
   seita_tmp = sort(seita_tmp);
   seita = (-seita_tmp + seita);
end
phi_e = seita/num;  
disp('ESPRIT算法信号源估计方向为：');
disp(phi_e);

end

