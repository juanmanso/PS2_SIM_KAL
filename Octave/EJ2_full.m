config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Estimacion a partir de mediciones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;			% Se considera solo x e y
tipos_variables = 3;		% Posicion, Velocidad, Aceleracion
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;

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

Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,dim)*var_xia]); %Solo para x e y



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cov_p = [1 1]*100^6;
cov_v = [1 1]*100;
cov_a = [1 1]*10;

x0 = [40 -200 0 0 0 0]';
P0_0 = diag([cov_p, cov_v, cov_a]);

%%%%%
sigma_etap = 100;
sigma_etav = 10;
sigma_etaa = 1;

%%% Para hacer AWGN, randn(fila,col)*sigma_eta

Bk1 = eye(cant_estados);

C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
	 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
	 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];

M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; 
	randn(dim,cant_mediciones)*sigma_etav*bool_v;
       	randn(dim,cant_mediciones)*sigma_etaa*bool_a];

yk = C * [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]' + (C*M_eta);
yk = yk'; % Asi tiene la forma de Pos

R = diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);


%%% ALGORITMO %%%%
x = x0;
P = P0_0;
xk1_k1 = x;
Pk1_k1 = P;
g = yk(1,:)';

for i=1:cant_mediciones-1
	% Prediccion
	xk_k1 = Ad * xk1_k1;
	Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1.';
	gk = [innovaciones(yk(i,:),C,xk_k1)];

	% Correccion
	Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
	xk_k = xk_k1 + Kk*(gk);
	Pk_k = (eye(cant_estados) - Kk*C) * Pk_k1;
	
	% Actualizacion
	xk1_k1 = xk_k;
	Pk1_k1 = Pk_k;


	% Guardo
	g = [g gk];
	x = [x xk_k];
	P = [P; Pk_k];
end


% Grafico de medida, estimada, ruidosa
h=figure;
subplot(2,2,1);
hold on
grid
if bool_p
    plot(yk(:,1),yk(:,2),'color',myGreen)
end
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x(1,:),x(2,:),'--b','LineWidth',2)
title('Estimacion de la trayectoria');
if(EsMatlab == 1)
    if bool_p
        legend('Medicion','Real','Estimada','location','SouthEast');
    else
        legend('Real','Estimada','location','SouthEast');
    end
    xlabel('Posicion x');
    ylabel('Posicion y');
else
    if bool_p
        legend(['Medicion';'Real';'Estimada'],'location','SouthEast');
    else
        legend(['Real';'Estimada'],'location','SouthEast');
    end
    xlabel('Posicion $x$ [\si{\m}]');
    ylabel('Posicion $y$ [\si{\m}]');
end

% Grafico del estado posicion en funcion del tiempo
%figure
subplot(2,2,2);
hold on
grid
plot(Pos(:,1),'LineWidth',2)
plot(Pos(:,2),'LineWidth',2)
plot(x(1,:),'--','LineWidth',2)
plot(x(2,:),'--','color',myGreen,'LineWidth',2)
ylabel('Posicion');
xlabel('Tiempo');
title('Estados de posicion');
legend('Real x','Real y','Estimada x','Estimada y','location','SouthEast');


% Grafico del estado velocidad en funcion del tiempo
% figure
subplot(2,2,3);
hold on
grid
if bool_v
    plot(yk(:,1),'--');
    plot(yk(:,2),'--');
end
plot(Vel(:,1),'LineWidth',2)
plot(Vel(:,2),'LineWidth',2)
plot(x(3,:),'--','LineWidth',2)
plot(x(4,:),'--','color',myGreen,'LineWidth',2)
ylabel('Velocidad');
xlabel('Tiempo');
title('Estados de velocidad');
if bool_v
    legend('Medida x','Medida y','Real x','Real y','Estimada x','Estimada y','location','SouthEast');
else
    legend('Real x','Real y','Estimada x','Estimada y','location','SouthEast');
end


% Grafico del estado aceleracion en funcion del tiempo
% figure
subplot(2,2,4);
hold on
grid
if bool_a
    plot(yk(:,1),'--');
    plot(yk(:,2),'--');
end
plot(Acel(:,1),'LineWidth',2)
plot(Acel(:,2),'LineWidth',2)
plot(x(5,:),'--','LineWidth',2)
plot(x(6,:),'--','color',myGreen,'LineWidth',2)
ylabel('Aceleracion');
xlabel('Tiempo');
title('Estados de aceleracion');
if bool_a
    legend('Medida x','Medida y','Real x','Real y','Estimada x','Estimada y','location','SouthEast');
else
    legend('Real x','Real y','Estimada x','Estimada y','location','SouthEast');
end

% Observabilidad
Obs = obsv(Ad,C);
rango_obs = rank(Obs);
estados_no_observables = cant_estados - rango_obs



