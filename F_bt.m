clear all
clc
N = 256;
Nfft = 256;
n = 0 : 1 : (N-1);
fai1 = 2*pi*rand(1,N);
fai2 = 2*pi*rand(1,N);
fai3 = 2*pi*rand(1,N);
s1 = sqrt(2000)*cos(2*pi*0.1*n+fai1);
s2 = sqrt(2000)*cos(2*pi*0.25*n+fai2);
s3 = sqrt((10)^(2.7)*2)*cos(2*pi*0.27*n+fai3);
v = rand(1,N);
u = s1 + s2 + s3 + v;

r_u = xcorr(u , 'biased'); 
S_bt = fft(r_u , Nfft);
S_bt = abs(S_bt);
S_bt=10*log10(S_bt);
index = n/N;
subplot(1,2,1);
plot(u);
subplot(1,2,2);
plot(index ,S_bt);
