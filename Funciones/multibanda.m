
% La siguiente función se encarga de generar un filtro
% con M dado, tal que entre las frecuencias de corte f
% (éste puede ser vector) amplifique con una ganancia a

%% a y f son vectores.
%% M es un escalar.
%% win es un "puntero a función" de ventanas.

% Si arrancase con un pasabanda (es decir, 
% suprimiendo las bajas frecuencias), el vector
% de frecuencias debe comenzar con 0.

%% Si win == @kaiser, la variable M contendrá M = [A, dw]

function [h, M] = multibanda(a,wc,M,win)

%	w = 2*pi*f;
	wc = [wc 1].*pi;
	a = [a];
	
	if(win == @kaiser)
		A = M(1); dw = M(2);
		[wd, M] = kaiser(A, dw);
	else
		wd = window(win,M+1)';
	end

	n = 0:M;

	% Inicializo hd
	hd = (0:M).*0;

	for i = 1:(length(wc))
		hd = hd + a(i).*sinc(wc(i)/pi*(n-M/2))*wc(i)/pi;
		if(i!=1)
			hd = hd - a(i).*sinc(wc(i-1)/pi*(n-M/2))*wc(i-1)/pi;
		end
	end


	h = hd .* wd;

end
