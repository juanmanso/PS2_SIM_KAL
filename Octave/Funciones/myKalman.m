% ALGORITMO DE KALMAN %

function [x,P,g] = myKalman(xk1_k1, Pk1_k1, Ak1, Bk1, Ck, Qk1, Rk, yk)

	% Predicción
	xk_k1 = Ak1 * xk1_k1;
	Pk_k1 =	Ak1 * Pk1_k1 * Ak1' + Bk1 * Qk1 * Bk1';
	gk = innovaciones(yk, Ck, xk_k1);

	% Corrección
	Kk = Pk_k1 * Ck'*(Rk + Ck*Pk_k1*Ck')^-1;
	xk_k = xk_k1 + Kk*(gk);
	Pk_k = (eye(length(Kk)) - Kk*Ck) * Pk_k1;
	
% PARA HACER SIMETRICA P, ALEJANDOSE DEL VALOR VERDADERO	
%		% Para hacerlo simétrico
%		for i=1:2
%			Pk_k=(Pk_k+Pk_k')/2;
%		end

	% Actualización
%	xk1_k1 = xk_k;
%	Pk1_k1 = Pk_k;

	g = gk;
	x = xk_k;
	P = Pk_k;
end

