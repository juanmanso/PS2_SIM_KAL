config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Variaci�n de R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;			% Se considera s�lo x e y
tipos_variables = 3;		% Posici�n, Velocidad, Aceleraci�n
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;

% Variables de configuraci�n
factor_R = 1;
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

Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,dim)*var_xia]); %S�lo para x e y



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cov_p = [1 1]*100^6;
cov_v = [1 1]*100;
cov_a = [1 1]*10;

x0 = [40 -200 0 0 0 0]';
P0_0 = diag([cov_p, cov_v, cov_a]);

%% a)
%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
sigma_etap = 100;
sigma_etav = 10;
sigma_etaa = 1;

%%% Para hacer AWGN, randn(fila,col)*sigma_etap

Bk1 = eye(cant_estados);

C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
	 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
	 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];

M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; 
	randn(dim,cant_mediciones)*sigma_etav*bool_v;
       	randn(dim,cant_mediciones)*sigma_etaa*bool_a];

yk = C * [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]' + (C*M_eta);
yk = yk'; % As� tiene la forma de Pos

R = factor_R*diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);


%%% ALGORITMO %%%%
x = x0;
P = P0_0;
xk1_k1 = x;
Pk1_k1 = P;
g = yk(1,:)';

for i=1:cant_mediciones-1
	% Predicci�n
	xk_k1 = Ad * xk1_k1;
	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
	gk = [innovaciones(yk(i,:),C,xk_k1)];

	% Correcci�n
	Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
	xk_k = xk_k1 + Kk*(gk);
	Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
	
	% Actualizaci�n
	xk1_k1 = xk_k;
	Pk1_k1 = Pk_k;


	% Guardo
	g = [g gk];
	x = [x xk_k];
	P = [P; Pk_k];
end


% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x(1,:),x(2,:),'LineWidth',3)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)
title('Estimaci�n');
if(EsMatlab == 1)
    legend('Estimada','Medida','Ruidosa');
else
    legend(['Estimada';'Medida';'Ruidosa']);
end
xlabel = 'Tiempo [s]';
ylabel = 'Posici�n [m]';

% Grafico del estado posici�n en funci�n del tiempo
figure
hold on
grid
plot(x(1,:),'LineWidth',2)
plot(x(2,:),'color',myGreen,'LineWidth',2)
title('Estados de posici�n');


% Grafico del estado velocidad en funci�n del tiempo
figure
hold on
grid
plot(x(3,:),'LineWidth',2)
plot(x(4,:),'color',myGreen,'LineWidth',2)
title('Estados de velocidad');


% Grafico del estado aceleraci�n en funci�n del tiempo
figure
hold on
grid
plot(x(5,:),'LineWidth',2)
plot(x(6,:),'color',myGreen,'LineWidth',2)
title('Estados de aceleración');


% Gr�fico de correlaci�n de innovaciones (debe ser ruido blanco)
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
estados_no_observables = cant_estados - rango_obs


