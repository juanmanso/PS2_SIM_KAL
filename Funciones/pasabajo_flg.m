

%% Filtro pasabajos

function [hd] = pasabajo_flg(wc,M)

	wc = wc*pi;  % frecuencia de corte del filtro
	n = 0:M;      % tiempo discreto

	% Respuesta temporal del pasabajos con retardo de M/2.
	hd = sinc(wc/pi*(n-M/2))*wc/pi;

endfunction

