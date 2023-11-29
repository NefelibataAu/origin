function [X] = Signal_Generator(k, z, phi, SNRs, N)
% SIGNAL_GENERATOR : 产生复高斯随机信号，输入来波方向，对应的信噪比以及快拍数，输出信号
% 输入：k:波数 z:阵元分布
% 输入：phi:来波方向 SNRs:信噪比 N:快拍数
% 输出：X:接收信号

%%

K = length(phi);
M = length(z);

X = zeros(M, N);

for i = 1 : K
    S = randn(1, N) + 1j * randn(1, N);             % 随机复高斯信号
    A = exp( -1j*k*z*sin(phi(i)) );                 % 导向向量
    X_temp = A * S;                                 % 单个源接收信号
    X_temp = awgn(X_temp, SNRs(i), 'measured');     % 添加噪声
    X = X + X_temp;                                 % 总的接收信号
end

end

