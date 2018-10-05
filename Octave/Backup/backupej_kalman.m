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

Qd = diag([var_xip,var_xip,var_xiv,var_xiv,var_xia,var_xia]); %Sólo para x e y



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = [40 -200 0 0 0 0]';
P0_0 = diag([100^2 100^2, 1 1, 0.1 0.1]);

%% a)
%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
sigma_etap = 60;
sigma_etav = 2;
sigma_etaa = 0.1;
cant_mediciones = length(Pos);

%%% Para hacer AWGN, randn(fila,col)*sigma_etap

C = [eye(dim) zeros(dim) zeros(dim)];

yk = C * [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]' + randn(dim,cant_mediciones)*sigma_etap;
yk = yk'; % Así tiene la forma de Pos

R = eye(dim)*sigma_etap^2;


%%% ALGORITMO %%%%
xk1_k1 = x0;
Pk1_k1 = P0_0;
Bk1 = eye(dim*3);
i=1;
x = x0;
P = P0_0;

for i=1:cant_mediciones-1
	% Predicción
	xk_k1 = Ad * xk1_k1;
	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';

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
	x = [x xk_k];
	P = [P; Pk_k];
end


% Grafico de medida, estimada, ruidosa
figure
hold on
plot(x(1,:),x(2,:),'LineWidth',3)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)

% Observabilidad
cant_estados = 6;
Obs = obsv(Ad,C);
rango_obs = rank(Obs);
estados_no_observables = cant_estados - rango_obs;

%%%%%%%%%%%%%%%%%%%%%%%
%%%% Ejercicio 3
%%%%%%%%%%%%%%%%%%%%%%
x01 = [40 -200 0 0 0 0]';
x02 = [200 -3000 0 0 0 0]';
x03 = x01;
x04 = x02;
x0_vec = [x01 x02 x03 x04];

P0_01 = diag([100^4 100^4, 10^2 10^2, 10 10]);
P0_02 = P0_01;
P0_03 = diag([0.1 0.1, 1e-5 1e-5, 1e-7 1e-7]);
P0_04 = P0_03;

xk1_k1 = x01;
Pk1_k1 = P0_01;
x = x01;
for i = 1:cant_mediciones-1
	[xk1_k1, Pk1_k1] = myKalman(xk1_k1,Pk1_k1, Ad, Bk1, C, Qd, R, yk(i,:));
	x1 = [x xk1_k1];
end

xk1_k1 = x02;
Pk1_k1 = P0_02;
x = x02;
for i = 1:cant_mediciones-1
	[xk1_k1, Pk1_k1] = myKalman(xk1_k1,Pk1_k1, Ad, Bk1, C, Qd, R, yk(i,:));
	x2 = [x xk1_k1];
end
	
xk1_k1 = x03;
Pk1_k1 = P0_03;
x = x03;
for i = 1:cant_mediciones-1
	[xk1_k1, Pk1_k1] = myKalman(xk1_k1,Pk1_k1, Ad, Bk1, C, Qd, R, yk(i,:));
	x3 = [x xk1_k1];
end

xk1_k1 = x04;
Pk1_k1 = P0_04;
x = x04;
for i = 1:cant_mediciones-1
	[xk1_k1, Pk1_k1] = myKalman(xk1_k1,Pk1_k1, Ad, Bk1, C, Qd, R, yk(i,:));
	x4 = [x xk1_k1];
end

% Grafico de medida, estimada, ruidosa
figure
hold on
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x(1,:),x(1,:),'color',[1 0 0]);
plot(x1(1,:),x1(1,:),'color',[0.8 0 0]);
plot(x2(1,:),x2(1,:),'color',[0.6 0 0]);
plot(x3(1,:),x3(1,:),'color',[0.4 0 0]);
plot(x4(1,:),x4(1,:),'color',[0.2 0 0]);



