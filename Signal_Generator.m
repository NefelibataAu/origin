function [X] = Signal_Generator(k, z, phi, SNRs, N)
% SIGNAL_GENERATOR : 产生窄带不相干随机信号，输入来波方向，对应的信噪比以及快拍数，输出信号
% 输入：k:波数 z:阵元分布
% 输入：phi:来波方向 SNRs:信噪比 N:快拍数
% 输出：X

%%

K = length(phi);
M = length(z);

X = zeros(M, N);

for i = 1 : K
    S = randn(1, N) + 1j * randn(1, N);
    A = exp( -1j*k*z*sin(phi(i)) );
    X_temp = awgn(A * S, SNRs(i), 'measured');
    X = X + X_temp;
end

end

