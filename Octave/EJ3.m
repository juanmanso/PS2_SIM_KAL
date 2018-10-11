config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Estimación con distintos valores iniciales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;			% Se considera sólo x e y
tipos_variables = 3;		% Posición, Velocidad, Aceleración
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

Qd = diag([ones(1,dim)*var_xip, ones(1,dim)*var_xiv,ones(1,dim)*var_xia]); %SÃ³lo para x e y



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

% C = [eye(dim)*bool_p eye(dim)*bool_v eye(dim)*bool_a];	% Obsoleto
C =	[eye(dim*bool_p) zeros(dim*bool_p) zeros(dim*bool_p);
	 zeros(dim*bool_v) eye(dim*bool_v) zeros(dim*bool_v);
	 zeros(dim*bool_a) zeros(dim*bool_a) eye(dim*bool_a)];

M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; 
	randn(dim,cant_mediciones)*sigma_etav*bool_v;
       	randn(dim,cant_mediciones)*sigma_etaa*bool_a];

yk = C * [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)]' + (C*M_eta);
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


% Nuevos datos
x01 = [40 -200 0 0 0 0]';
x02 = [200 -3000 0 0 0 0]';
x03 = x01;
x04 = x02;
x0_vec = [x01 x02 x03 x04];

diag1_p = 100*[cov_p, cov_v, cov_a];
diag2_p = 0.01*[cov_p, cov_v, cov_a];
diag_vec = [diag1_p; diag1_p; diag2_p; diag2_p];


for l = 1:4
	xa = x0_vec(:,l);
	Pa = diag(diag_vec(l,:));
	xk1_k1 = xa;
	Pk1_k1 = Pa;

	for i = 1:cant_mediciones-1
		% Predicción
		xk_k1 = Ad * xk1_k1;
		Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
	
		% Corrección
		Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
		xk_k = xk_k1 + Kk*(yk(i,:)' - C*xk_k1);
		Pk_k = (eye(dim*3) - Kk*C) * Pk_k1;
		
	
		% Actualización
		xk1_k1 = xk_k;
		Pk1_k1 = Pk_k;
	
	
		% Guardo
		xa = [xa xk_k];
		Pa = [Pa; Pk_k];
	end

	if(l==1)
		x1=xa; P1=Pa;
	elseif(l==2)
		x2=xa; P2=Pa;
	elseif(l==3)
		x3=xa; P3=Pa;
	else
		x4=xa; P4=Pa;
	end
	
end


% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x(1,:),x(2,:),'LineWidth',3)
%	plot(x1(1,:),x1(2,:),'-','LineWidth',3,'color',[0 0.8 1])
%	plot(x2(1,:),x2(2,:),'-','LineWidth',3,'color',[0.8 0.8 1])
%	plot(x3(1,:),x3(2,:),'-','LineWidth',3,'color',[0.4 0.4 0.4])
%	plot(x4(1,:),x4(2,:),'-','LineWidth',3,'color',[0.8 0.2 1])
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)
title('Estimación');
if(EsMatlab == 1)
    legend('Estimada','Medida','Ruidosa','location','SouthEast');
else
    legend(['Estimada';'Medida';'Ruidosa'],'location','SouthEast');
end
xlabel=('Posición $x$ [\si{\m}]');
ylabel=('Posición $y$ [\si{\m}]');


% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x1(1,:),x1(2,:),'LineWidth',3)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)
title('Estimación');
if(EsMatlab == 1)
    legend('Estimada','Medida','Ruidosa','location','SouthEast');
else
    legend(['Estimada';'Medida';'Ruidosa'],'location','SouthEast');
end
xlabel=('Posición $x$ [\si{\m}]');
ylabel=('Posición $y$ [\si{\m}]');


% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x2(1,:),x2(2,:),'LineWidth',3)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)
title('Estimación');
if(EsMatlab == 1)
    legend('Estimada','Medida','Ruidosa','location','SouthEast');
else
    legend(['Estimada';'Medida';'Ruidosa'],'location','SouthEast');
end
xlabel=('Posición $x$ [\si{\m}]');
ylabel=('Posición $y$ [\si{\m}]');

% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x3(1,:),x3(2,:),'LineWidth',3)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)
title('Estimación');
if(EsMatlab == 1)
    legend('Estimada','Medida','Ruidosa','location','SouthEast');
else
    legend(['Estimada';'Medida';'Ruidosa'],'location','SouthEast');
end
xlabel=('Posición $x$ [\si{\m}]');
ylabel=('Posición $y$ [\si{\m}]');

% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(x4(1,:),x4(2,:),'LineWidth',3)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(yk(:,1),yk(:,2),'color',myGreen)
title('Estimación');
if(EsMatlab == 1)
    legend('Estimada','Medida','Ruidosa','location','SouthEast');
else
    legend(['Estimada';'Medida';'Ruidosa'],'location','SouthEast');
end
xlabel=('Posición $x$ [\si{\m}]');
ylabel=('Posición $y$ [\si{\m}]');


%%%%%%%%%% CÃ³mo grafico el error?!?! %%%%%%%%%%%%%%%%%%
