config_m;
%%%%%%%%%%%%%%
% EJ KALMAN
%%%%%%%%%%%%%%

% Algortimo de Filtro de Kalman
% Conozco:	- dinámica (A,B,C,D)
%		- Parámetros Estadísticos \psi, \eta
%		- Condiciones iniciales: \hat{X}_{0|0} = \mathbb{E}(X_o), P_{0|0} = cov (X_{0}_\hat{X}_{0|0})

% Predicción:	\hat{X}_{k|k-1} = A_{k-1} \cdot \hat{X}_{k-1|k-1} 
%		P_{k|k-1} = A_{k-1} P_{k-1|k-1} A^{*}_{k-1} + B_{k-1} Q_{k-1} B^{*}_{k-1}
% Corrección:	K_k = P_{k|k-1}

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;

bool_p = 1;
bool_v = 0;
bool_a = 0;

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
T = 1;

% Variable de estado X = [P;V;A]
I = eye(dim);
Ad =	[I	I.*T	(T.^2)/2.*I;
	 I*0	I	T.*I;
	 I*0	I*0	I;];

Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,dim)*var_xia]); %Sólo para x e y



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJERCICIO 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [40 -200 0 0 0 0]';
pb = [30 45]';
vb = [0 0]';
ab = [0 0]';
P0_0 = diag([100^2 100^2, 1 1, 0.1 0.1]);

%% a)
%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
sigma_etap = 60;
sigma_etav = 2;
sigma_etaa = 0.1;
cant_mediciones = length(Pos);

%%% Para hacer AWGN, randn(fila,col)*sigma_etap

C = [eye(dim)*bool_p eye(dim)*bool_v eye(dim)*bool_a];
b = [ones(cant_mediciones,dim).*pb' ones(cant_mediciones,dim).*vb' ones(cant_mediciones,dim).*ab'];

yk = C * ([Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)])' + randn(dim,cant_mediciones)*sigma_etap;
yk = yk'; % Así tiene la forma de Pos

R = eye(dim)*sigma_etap^2;


%%% ALGORITMO %%%%
xk1_k1 = x0;
Pk1_k1 = P0_0;
Bk1 = eye(dim*3);
i=1;
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

% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x(1,:),x(2,:),'LineWidth',3)
plot(xb(1,:),xb(2,:),'LineWidth',3, 'color', myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color','k')
title('Estimación');
legend(['Estimada sin corregir'; 'Estimada corregida';'Medida';'Ruidosa']);
xlabel = 'Tiempo [s]';
ylabel = 'Posición [m]';


% Grafico del estado posición en función del tiempo
figure
hold on
grid
plot(x(1,:),'LineWidth',2)
plot(xb(1,:),'color',myGreen,'LineWidth',2)
legend(['Sin corregir';'Corregida'])
title('Estados de posición en x');

% Grafico del estado posición en función del tiempo
figure
hold on
grid
plot(x(2,:),'LineWidth',2)
plot(xb(2,:),'color',myGreen,'LineWidth',2)
legend(['Sin corregir';'Corregida'])
title('Estados de posición en y');

return;


% Grafico del estado velocidad en función del tiempo
figure
hold on
grid
plot(x(3,:),'LineWidth',2)
plot(x(4,:),'color',myGreen,'LineWidth',2)
title('Estados de velocidad');


% Grafico del estado aceleración en función del tiempo
figure
hold on
grid
plot(x(5,:),'LineWidth',2)
plot(x(6,:),'color',myGreen,'LineWidth',2)
title('Estados de aceleración');


% Gráfico de Ruido Blanco
figure
plot(g(1,:),g(2,:));
mean(g') 
sqrt(var(g')) 

% Observabilidad
cant_estados = 3*dim;
Obs = obsv(Ad,C);
rango_obs = rank(Obs);
estados_no_observables = cant_estados - rango_obs


