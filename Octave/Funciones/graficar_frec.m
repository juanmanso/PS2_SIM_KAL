

function [] = graficar_frec(hd)

nfft = 1*1024;  % tamaño de la FFT 
w = (0:nfft/2-1)/nfft*pi*2;

H = fft(hd,nfft); H = H(1:end/2);

plot(w/pi,abs(H));

grid minor

endfunction
