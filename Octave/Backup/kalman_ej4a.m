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

%datos_str = load('datos.mat');
%
%Acel = datos_str.Acel;
%Tiempo = datos_str.tiempo;
%Pos = datos_str.Pos;
%Vel = datos_str.Vel;

load('datos.mat');

dim = 2;			% Se considera sólo x e y
tipos_variables = 3;		% Posición, Velocidad, Aceleración
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;

bool_p = 1;
bool_v = 0;
bool_a = 0;

bool_pb = 1;
bool_vb = 0;
bool_ab = 0;
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
%T = Tiempo(2:end)-Tiempo(1:end-1);	% Si A fuese dinámica en función de k
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

cov_p = [1 1]*100^2;
cov_v = [1 1]*1;
cov_a = [1 1]*0.1;

%x0 = [40 -200 0 0 0 0]';
x0 = [10 -50 0 0 0 0]';
P0_0 = diag([cov_p, cov_v, cov_a]);

%% a)
%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
sigma_etap = 60;
sigma_etav = 2;
sigma_etaa = 0.1;

%%% Para hacer AWGN, randn(fila,col)*sigma_etap

Bk1 = eye(cant_estados);

C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
	 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
	 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];

M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; randn(dim,cant_mediciones)*sigma_etav*bool_v; randn(dim,cant_mediciones)*sigma_etaa*bool_a];


yk = C * ([Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)])' + (C*M_eta);
yk = yk'; % Así tiene la forma de Pos

R = diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);


%%% ALGORITMO %%%%
x = x0;
P = P0_0;
xk1_k1 = x;
Pk1_k1 = P;
g = yk(1,:)';

for i=1:cant_mediciones-1
	% Predicción
	xk_k1 = Ad * xk1_k1;
	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
	gk = [innovaciones(yk(i,:),C,xk_k1)];

	% Corrección
	Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
	xk_k = xk_k1 + Kk*(gk);
	Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
	
	% Actualización
	xk1_k1 = xk_k;
	Pk1_k1 = Pk_k;

	% Guardo
	g = [g gk];
	x = [x xk_k];
	P = [P; Pk_k];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Con corrección de sesgo %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supongo que desconozco un sesgo a la vez
var_xib = 0;
b0 = [0 0]';
cov_b = cov_p*bool_pb + cov_v*bool_vb + cov_a*bool_ab;

% Redefino A y C suponiendo var de estado z = [x; b]
Ad_b = [Ad zeros(cant_estados,dim); zeros(dim,cant_estados) eye(dim)];
B_b = diag([ones(1,dim), ones(1,dim), ones(1,dim), zeros(1,dim)]);	% Zeros porque la dinámica del sesgo es cte
Qd_b = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv, ones(1,dim)*var_xia, ones(1,dim)*var_xib]); 
C_b = [C [eye(dim*bool_p)*bool_pb; eye(dim*bool_v)*bool_vb; eye(dim*bool_a)*bool_ab]]; % Concateno columnas de 0's excepto identidad en donde corresponda

tipos_variables = 4;
cant_estados = tipos_variables * dim;

% Problema con el eta
Sesgo = zeros(cant_mediciones,dim);

yk = C_b * ([Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim) Sesgo])' + C_b * [(M_eta); zeros(dim,cant_mediciones)];
yk = yk';

% ALGORITMO
xb = [x0; b0];
Pb = diag([cov_p, cov_v, cov_a, cov_b]);
xk1_k1 = xb;
Pk1_k1 = Pb;
gb = yk(1,:)';


for i=1:cant_mediciones-1
	% Predicción
	xk_k1 = Ad_b * xk1_k1;
	Pk_k1 =	Ad_b * Pk1_k1 * Ad_b' + B_b * Qd_b * B_b';
	gk = [innovaciones(yk(i,:),C_b,xk_k1)];

	% Corrección
	Kk = Pk_k1 * C_b'*(R + C_b*Pk_k1*C_b')^-1;
	xk_k = xk_k1 + Kk*(gk);
	Pk_k = (eye(cant_estados) - Kk*C_b) * Pk_k1;
	
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
%plot(yk(:,1),yk(:,2),'color','k')
title('Estimación');
legend(['Estimada sin corregir'; 'Estimada corregida';'Medida';'Ruidosa'],'location','SouthEast');
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

% Gráfico de correlación de innovaciones (debe ser ruido blanco)
covx_g = xcorr(g(1,:)');
covy_g = xcorr(g(2,:)');

figure
plot(covx_g)
grid
title('Covarianza innovaciones x')

figure
plot(covy_g)
grid
title('Covarianza innovaciones y')


% Observabilidad
Obs = obsv(Ad,C);
rango_obs = rank(Obs);
estados_no_observables = (cant_estados) - rango_obs


