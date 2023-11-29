function [] = draw_dB(P_dB)
% DRAW_DB 给定空间谱，进行作图
% 输入：P_dB:空间谱
P_dB = P_dB';
[Num, Len] = size(P_dB);

figure;

seita = linspace(-90,90,Len)*pi/180;
for i = 1:Num
    plot(seita*180/pi, P_dB(i,:), 'k', 'Linewidth', 2);
    xlabel('空间角度 (deg)');
    ylabel('空间谱 (dB)');
    grid on;
    hold on;
end

end

