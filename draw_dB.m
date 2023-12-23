function [] = draw_dB(P_dB)
% DRAW_DB 给定空间谱，进行作图
% 输入：P_dB:空间谱

[Num, Len] = size(P_dB);

figure;

seita = linspace(-90,90,Len);
for i = 1:Num
    plot(seita, P_dB(i,:));
    hold on;
end
xlabel('空间角度 (deg)');
ylabel('空间谱 (dB)');
grid on;

end

