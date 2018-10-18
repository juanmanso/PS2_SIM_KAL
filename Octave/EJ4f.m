config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Estimación Sesgo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datos_str = load('datos.mat');

acel = datos_str.Acel;
tiempo = datos_str.tiempo;
pos = datos_str.Pos;
Vel = datos_str.Vel;

load('datos.mat');

dim = 2;			% Se considera sólo x e y
tipos_variables = 3;		% Posición, Velocidad, Aceleración
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;

% Selección de medición (se pueden múltiples opciones)
bool_p = 1;
bool_v = 1;
bool_a = 1;

% Selección de sesgo (sólo 1)
bool_pb = 0;
bool_vb = 0;
bool_ab = 1;

% Selección de impresión de imágenes
bool_print = 1;

% Supongo que desconozco un sesgo a la vez
b0_p = [300 200]';
b0_v = [10 20]';
b0_a = [2 1]';
var_xib = 0;

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

% Con esta estructura las matrices tienen dimensiones dinÃ¡micas
C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
	 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
	 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];

% A partir de qué medición entra, se le agrega el ruido correspondiente a ella.
M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p+repmat(b0_p,1,cant_mediciones)*bool_pb; randn(dim,cant_mediciones)*sigma_etav*bool_v+repmat(b0_v,1,cant_mediciones)*bool_vb; randn(dim,cant_mediciones)*sigma_etaa*bool_a+repmat(b0_a,1,cant_mediciones)*bool_ab];


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

b0 = b0_p*bool_pb + b0_v*bool_vb + b0_a*bool_ab;
% cov_b = cov_p*bool_pb + cov_v*bool_vb + cov_a*bool_ab;
cov_b=[200^2 200^2];

% Redefino A y C suponiendo var de estado z = [x; b]
Ad_b = [Ad zeros(cant_estados,dim); zeros(dim,cant_estados) eye(dim)];
B_b = diag([ones(1,dim), ones(1,dim), ones(1,dim), zeros(1,dim)]);	% Zeros porque la dinámica del sesgo es cte
Qd_b = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv, ones(1,dim)*var_xia, ones(1,dim)*var_xib]); 
C_b = [C [eye(dim*bool_p)*bool_pb; eye(dim*bool_v)*bool_vb; eye(dim*bool_a)*bool_ab]]; % Concateno columnas de 0's excepto identidad en donde corresponda
R = diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);

tipos_variables = 4;
cant_estados = tipos_variables * dim;

% Problema con el eta
Sesgo = zeros(cant_mediciones,dim);

yk = C_b * ([Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim) Sesgo])' + C_b * [(M_eta); zeros(dim,cant_mediciones)];
yk = yk';

% ALGORITMO
% xb = [x0; b0];
xb = [x0; 0;0];
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
h1=figure;
subplot(221)
hold on
grid
plot(yk(:,1),yk(:,2))
plot(x(1,:),x(2,:),'m','LineWidth',2)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(xb(1,:),xb(2,:),'--b','LineWidth',2)
title('Estimacion');
if(EsMatlab == 1)
    legend('Medido','Estimación sin estimar sesgo','Estado','Estimada estimando sesgo','location','SouthEast');
else 
    legend(['Medido';'Estimación sin estimar sesgo';'Estado';'Estimada estimando sesgo'],'location','SouthEast');
end
xlabel = 'Tiempo [s]';
ylabel = 'Posición [m]';


% Grafico del estado posición en función del tiempo
% figure
subplot(222)
hold on
grid
plot(x(1,:),'LineWidth',2)
plot(xb(1,:),'--','LineWidth',2)
plot(Pos(:,1));
plot(x(2,:),'LineWidth',2)
plot(xb(2,:),'--','LineWidth',2)
plot(Pos(:,2));
if(EsMatlab == 1)
    legend('px s/sesgo','px c/sesgo','Estado px','py s/sesgo','py c/sesgo','Estado py','location','SouthWest')
else 
    legend(['px s/sesgo';'px c/sesgo';'Estado px';'py s/sesgo';'py c/sesgo';'Estado py'],'location','SouthWest')
end
title('Estados de posición');

% Grafico del estado velocidad en función del tiempo
% figure
subplot(223)
hold on
grid
plot(x(3,:),'LineWidth',2)
plot(xb(3,:),'--','LineWidth',2)
plot(Vel(:,1),'b','LineWidth',1);
plot(x(4,:),'LineWidth',2)
plot(xb(4,:),'--','color',myGreen,'LineWidth',2)
plot(Vel(:,2),'r','LineWidth',1);

if(EsMatlab == 1)
    legend('Vx s/sesgo','Vx c/sesgo','Estado Vx','Vy s/sesgo','Vy c/sesgo','Estado Vy','location','SouthEast')
else 
    legend(['Vx s/sesgo';'Vx c/sesgo';'Estado Vx';'Vy s/sesgo';'Vy c/sesgo';'Estado Vy'],'location','SouthEast')
end
title('Estados de velocidad');

% Grafico del estado aceleración en función del tiempo
% figure
subplot(224)
hold on
grid
plot(x(5,:),'LineWidth',2)
plot(xb(5,:),'--','LineWidth',2)
plot(Acel(:,1),'b','LineWidth',1);
plot(x(6,:),'LineWidth',2)
plot(xb(6,:),'--','color',myGreen,'LineWidth',2)
plot(Acel(:,2),'r','LineWidth',1);

if(EsMatlab == 1)
    legend('Ax s/sesgo','Ax c/sesgo','Estado Ax','Ay s/sesgo','Ay c/sesgo','Estado Ay','location','SouthEast')
else
    legend(['Ax s/sesgo';'Ax c/sesgo';'Estado Ax';'Ay s/sesgo';'Ay c/sesgo';'Estado Ay'],'location','SouthEast')
end
title('Estados de aceleración');

h1.Position=[0 0 1200 700];
h1.PaperUnits='points';
h1.PaperSize=[1200 700];
if bool_print
    print('../Informe/Figuras/graf_ej4f','-dpdf','-bestfit');
end

% Grafico de estimación de sesgo
h2=figure;
hold on
grid
plot(xb(7,:),'LineWidth',3)
plot(xb(8,:),'LineWidth',3)
title('Estimacion de sesgo');
if(EsMatlab == 1)
    legend('Sesgo x','Sesgo y','location','SouthEast');
else 
    legend(['Sesgo x','Sesgo y'],'location','SouthEast');
end
xlabel = 'Tiempo [s]';
ylabel = 'Posición [m]';

wsize=[h2.Position(3) h2.Position(4)];
h2.PaperUnits='points';
h2.PaperSize=wsize;
if bool_print
    print('../Informe/Figuras/bias_ej4f','-dpdf','-bestfit');
end


% Gráfico de correlación de innovaciones (debe ser ruido blanco)
covx_g = xcorr(g(1,:)');
covy_g = xcorr(g(2,:)');

h3=figure;
subplot(211)
plot(covx_g)
grid
axis tight;
title('Covarianza innovaciones x')

% figure
subplot(212)
plot(covy_g)
grid
axis tight;
title('Covarianza innovaciones y')

wsize=[h3.Position(3) h3.Position(4)];
h3.PaperUnits='points';
h3.PaperSize=wsize;
if bool_print
    print('../Informe/Figuras/covinn_ej4f','-dpdf','-bestfit');
end

% Observabilidad
Obs = obsv(Ad_b,C_b);
rango_obs = rank(Obs);
estados_no_observables = (cant_estados) - rango_obs
