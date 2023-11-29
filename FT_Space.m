function [P_FT_dB] = FT_Space(k, z, X)
% FT_SPACE 空间傅里叶变换求DOA
% 输入：k:波数 z:阵元分布
% 输入：X:接收信号
% 输出：P_FT_dB:空间谱功率

seita = linspace(-pi/2, pi/2, 500);
a_seita = exp(-1j*k*z*sin(seita));
P_FT = a_seita' * X; 

maxdB = max( sum( abs(P_FT).^2, 2 ) );
P_FT_dB = 10 * log10( sum( (abs(P_FT).^2 )/maxdB, 2 ) ) ;

end

