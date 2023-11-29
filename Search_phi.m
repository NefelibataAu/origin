function [phi_e] = Search_phi(P_dB, K)
% Search_phi 通过信号空间谱（伪谱）求信号来向
% 输入： P_dB：空间谱 K:信号源个数
% 输出： phi：方向

len = length(P_dB);
seita = linspace(-90,90,len)*pi/180;

[P_peaks, P_peaks_idx] = findpeaks(P_dB);           % 提取峰值
[~, I] = sort(P_peaks, 'descend');                  % 峰值降序排序
P_peaks_idx = P_peaks_idx(I);                       % 提取前M个
P_peaks_idx = P_peaks_idx(1:K);
phi_e = seita(P_peaks_idx)*180/pi;                  % 估计方向

end

