config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Kalman Estacionario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;			% Se considera s�lo x e y
tipos_variables = 3;		% Posici�n, Velocidad, Aceleraci�n
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;

bool_p = 1;
bool_v = 1;
bool_a = 1;

bool_print = 0; % Flag para exportar imagen

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

x0 = [40 -1000 0 0 0 0]';
P0_0 = diag([cov_p, cov_v, cov_a]);

%% a)
%%%%% y_k = [I 0 0] [pk vk ak]' + ruido \eta
sigma_etap = 100;
sigma_etav = 10;
sigma_etaa = 1;

%%% Para hacer AWGN, randn(fila,col)*sigma_etap

Bk1 = eye(cant_estados);

% C = [eye(dim)*bool_p eye(dim)*bool_v eye(dim)*bool_a];	% Obsoleto
C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
	 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
	 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];

M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; 
	randn(dim,cant_mediciones)*sigma_etav*bool_v;
       	randn(dim,cant_mediciones)*sigma_etaa*bool_a];

yk = C * [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]' + (C*M_eta);
yk = yk'; % As� tiene la forma de Pos

R = diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);

%%% Funci�n DARE %%%

[P_dare, Polos_Lazo_Cerrado, K_dare] = dare(Ad,Bk1,Qd,R);

%%% ALGORITMO %%%%
x = x0;
P = P_dare;
xk1_k1 = x;
Pk1_k1 = P;
g = yk(1,:)';

for i=1:cant_mediciones-1
	% Predicci�n
	xk_k1 = Ad * xk1_k1;
% 	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
	gk = [innovaciones(yk(i,:),C,xk_k1)];

	% Correcci�n
% 	Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
	xk_k = xk_k1 + K_dare*(gk);
% 	Pk_k = (eye(cant_estados) - K_dare*C) * Pk_k1;
	
	% Actualizaci�n
	xk1_k1 = xk_k;
% 	Pk1_k1 = Pk_k;


	% Guardo
	g = [g gk];
	x = [x xk_k];
% 	P = [P; Pk_k];
end

% Grafico de medida, estimada, ruidosa
h=figure;
subplot(2,2,1);
hold on
grid
plot(yk(:,1),yk(:,2),'color',myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x(1,:),x(2,:),'--b','LineWidth',2)
title('Estimaci�n de la trayectoria');
if(EsMatlab == 1)
    legend('Medici�n','Real','Estimada','location','SouthEast');
    xlabel('Posici�n x');
    ylabel('Posici�n y');
else
    legend(['Medici�n';'Real';'Estimada'],'location','SouthEast');
    xlabel('Posici�n $x$ [\si{\m}]');
    ylabel('Posici�n $y$ [\si{\m}]');
end

% Grafico del estado posici�n en funci�n del tiempo
%figure
subplot(2,2,2);
hold on
grid
plot(Pos(:,1),'LineWidth',2)
plot(Pos(:,2),'LineWidth',2)
plot(x(1,:),'--','LineWidth',2)
plot(x(2,:),'--','color',myGreen,'LineWidth',2)
ylabel('Posici�n');
xlabel('Tiempo');
title('Estados de posici�n');
legend('Real x','Real y','Estimada x','Estimada y','location','SouthEast');


% Grafico del estado velocidad en funci�n del tiempo
% figure
subplot(2,2,3);
hold on
grid
plot(Vel(:,1),'LineWidth',2)
plot(Vel(:,2),'LineWidth',2)
plot(x(3,:),'--','LineWidth',2)
plot(x(4,:),'--','color',myGreen,'LineWidth',2)
plot(yk(:,3));
plot(yk(:,4));
ylabel('Velocidad');
xlabel('Tiempo');
title('Estados de velocidad');
legend('Real x','Real y','Estimada x','Estimada y','location','SouthEast');


% Grafico del estado aceleraci�n en funci�n del tiempo
% figure
subplot(2,2,4);
hold on
grid
plot(Acel(:,1),'LineWidth',2)
plot(Acel(:,2),'LineWidth',2)
plot(x(5,:),'--','LineWidth',2)
plot(x(6,:),'--','color',myGreen,'LineWidth',2)
plot(yk(:,5));
plot(yk(:,6));
ylabel('Aceleraci�n');
xlabel('Tiempo');
title('Estados de aceleraci�n');
legend('Real x','Real y','Estimada x','Estimada y','location','SouthEast');

h.Position=[0 0 1200 700];
h.PaperUnits='points';
h.PaperSize=[1200 700];
if bool_print
    print('../Informe/Figuras/graf_ej5','-dpdf','-bestfit');
end

% Gr�fico de correlaci�n de innovaciones (debe ser ruido blanco)
covx_g = xcorr(g(1,:)');
covy_g = xcorr(g(2,:)');

h2=figure;
subplot(211)
plot(covx_g)
grid
title('Covarianza innovaciones x')
axis tight;

% figure
subplot(212)
plot(covy_g)
grid
title('Covarianza innovaciones y')
axis tight;

h2.Position=[0 0 1200 700];
h2.PaperUnits='points';
h2.PaperSize=[1200 700];
if bool_print
    print('../Informe/Figuras/covinn_ej5','-dpdf','-bestfit');
end
