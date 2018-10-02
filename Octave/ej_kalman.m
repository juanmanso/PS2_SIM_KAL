config_m;
%%%%%%%%%%%%%%
% EJ KALMAN
%%%%%%%%%%%%%%

% Algortimo de Filtro de Kalman
% Conozco:	- dinámica (A,B,C,D)
%		- Parámetros Estadísticos \xi, \eta
%		- Condiciones iniciales: \hat{X}_{0|0} = \mathbb{E}(X_o), P_{0|0} = cov (X_{0}_\hat{X}_{0|0})

% Predicción:	\hat{X}_{k|k-1} = A_{k-1} \cdot \hat{X}_{k-1|k-1} 
%		P_{k|k-1} = A_{k-1} P_{k-1|k-1} A^{*}_{k-1} + B_{k-1} Q_{k-1} B^{*}_{k-1}
% Corrección:	K_k = P_{k|k-1}

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;			% Se considera sólo x e y
tipos_variables = 3;		% Posición, Velocidad, Aceleración
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;


%%%%%%%%%%%%%%
%%% 1a Defina las variables de estado
%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
%%% 1b Ad y Q_d
%%%%%%%%%%%%%%%

% Datos
var_xip = 3e-4;
var_xiv = 2e-3;
var_xia = 1e-2;

%%%
T = Tiempo(2:end)-Tiempo(1:end-1);
%T = 1;					% Suponiendo equiespaciado

% Variable de estado X = [P;V;A]
I = eye(dim);
Ad =	[I	I.*T	(T.^2)/2.*I;
	 I*0	I	T.*I;
	 I*0	I*0	I;];

% Covarianza del ruido de proceso
Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,dim)*var_xia]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJERCICIO 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bool_p = 1;	% Inciso a
bool_v = 0;	% Inciso b
bool_a = 0;	% Inciso c


x0 = [40 -200 0 0 0 0]';
P0_0 = diag([100^2 100^2, 1 1, 0.1 0.1]);

%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
sigma_etap = 60;
sigma_etav = 2;
sigma_etaa = 0.1;


%%% Para hacer AWGN, randn(fila,col)*sigma_etap

C = [eye(dim)*bool_p eye(dim)*bool_v eye(dim)*bool_a];

yk = C * ([Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)])' + randn(dim,cant_mediciones)*sigma_etap;
yk = yk'; % Así tiene la forma de Pos

R = eye(dim)*sigma_etap^2;
Bk1 = eye(dim*3);


%%% ALGORITMO %%%%
xk1_k1 = x0;
Pk1_k1 = P0_0;
x = x0;
P = P0_0;
g = yk(1,:)';

for i=1:cant_mediciones-1
	% Predicción
	xk_k1 = Ad * xk1_k1;
	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
	gk = [innovaciones(yk(i,:),C,xk_k1)];

	% Corrección
	Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
	xk_k = xk_k1 + Kk*(yk(i,:)' - C*xk_k1);
	Pk_k = (eye(dim*3) - Kk*C) * Pk_k1;
	
% PARA HACER SIMETRICA P, ALEJANDOSE DEL VALOR VERDADERO	
%		% Para hacerlo simétrico
%		for i=1:2
%			Pk_k=(Pk_k+Pk_k')/2;
%		end

	% Actualización
	xk1_k1 = xk_k;
	Pk1_k1 = Pk_k;


	% Guardo
	g = [g gk];
	x = [x xk_k];
	P = [P; Pk_k];
end

%% Con corrección de sesgo
xk1_k1 = x0;
Pk1_k1 = P0_0;
Bk1 = eye(dim*3);
i=1;
xb = x0;
Pb = P0_0;
gb = yk(1,:)';

yk = C * ([Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]-b)' + randn(dim,cant_mediciones)*sigma_etap;
yk = yk';

for i=1:cant_mediciones-1
	% Predicción
	xk_k1 = Ad * xk1_k1;
	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
	gk = [innovaciones(yk(i,:),C,xk_k1)];

	% Corrección
	Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
	xk_k = xk_k1 + Kk*(yk(i,:)' - C*xk_k1);
	Pk_k = (eye(dim*3) - Kk*C) * Pk_k1;
	
% PARA HACER SIMETRICA P, ALEJANDOSE DEL VALOR VERDADERO	
%		% Para hacerlo simétrico
%		for i=1:2
%			Pk_k=(Pk_k+Pk_k')/2;
%		end

	% Actualización
	xk1_k1 = xk_k;
	Pk1_k1 = Pk_k;


	% Guardo
	gb = [gb gk];
	xb = [xb xk_k];
	Pb = [Pb; Pk_k];
end

% Grafico del estado velocidad en función del tiempo
figure
hold on
grid
plot(x(3,:),'LineWidth',2)
plot(xb(3,:),'LineWidth',3, 'color', myGreen)
legend('Sin corrección','Corregida')
title('Estados de velocidad en x');

figure
hold on
grid
plot(x(4,:),'LineWidth',3, 'color', 'b')
plot(xb(4,:),'color',myGreen,'LineWidth',2)
legend('Sin corrección','Corregida')
title('Estados de velocidad en y');


