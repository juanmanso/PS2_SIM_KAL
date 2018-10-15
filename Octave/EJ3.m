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
    g = yk(1,:)';
    
	for i = 1:cant_mediciones-1
		% Predicción
		xk_k1 = Ad * xk1_k1;
		Pk_k1 =	Ad * Pk1_k1 * Ad' + Bk1 * Qd * Bk1';
        gk = [innovaciones(yk(i,:),C,xk_k1)];
	
		% Corrección
		Kk = Pk_k1 * C'*(R + C*Pk_k1*C')^-1;
		xk_k = xk_k1 + Kk*(yk(i,:)' - C*xk_k1);
		Pk_k = (eye(dim*3) - Kk*C) * Pk_k1;
		
	
		% Actualización
		xk1_k1 = xk_k;
		Pk1_k1 = Pk_k;
	
	
		% Guardo
        g = [g gk];
		xa = [xa xk_k];
		Pa = [Pa; Pk_k];
	end

	if(l==1)
		x1=xa; P1=Pa; g1=g;
	elseif(l==2)
		x2=xa; P2=Pa; g2=g;
	elseif(l==3)
		x3=xa; P3=Pa; g3=g;
	else
		x4=xa; P4=Pa; g4=g;
	end
	
end


% Grafico de medida, estimada, ruidosa
figure
hold on
grid
plot(yk(:,1),yk(:,2),'color',myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x(1,:),x(2,:),'--b','LineWidth',2)
%	plot(x1(1,:),x1(2,:),'-','LineWidth',3,'color',[0 0.8 1])
%	plot(x2(1,:),x2(2,:),'-','LineWidth',3,'color',[0.8 0.8 1])
%	plot(x3(1,:),x3(2,:),'-','LineWidth',3,'color',[0.4 0.4 0.4])
%	plot(x4(1,:),x4(2,:),'-','LineWidth',3,'color',[0.8 0.2 1])
title('Estimación');
if(EsMatlab == 1)
    legend('Medida','Real','Estimada','location','SouthEast');
    xlabel('Posición x');
    ylabel('Posición y');
else
    legend(['Medida';'Real','Estimada'],'location','SouthEast');
    xlabel('Posición $x$ [\si{\m}]');
    ylabel('Posición $y$ [\si{\m}]');
end

% Grafico de medida, estimada, ruidosa
h1=figure;
subplot(221)
hold on
grid
plot(yk(:,1),yk(:,2),'color',myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x1(1,:),x1(2,:),'--b','LineWidth',2)
title('Estimación');
if(EsMatlab == 1)
    legend('Medida','Real','Estimada','location','SouthEast');
    xlabel('Posición x');
    ylabel('Posición y');
else
    legend(['Medida';'Real','Estimada'],'location','SouthEast');
    xlabel('Posición $x$ [\si{\m}]');
    ylabel('Posición $y$ [\si{\m}]');
end

% figure;
subplot(222)
plot(x'-x1');
xlim([0 150]);
title('Error de estimación');
xlabel('Tiempo');
legend('e_p_x','e_p_y','e_v_x','e_v_y','e_a_x','e_a_y');

covx_g1 = xcorr(g1(1,:)');
covy_g1 = xcorr(g1(2,:)');

% figure
subplot(223)
plot(covx_g1)
grid
title('Covarianza innovaciones x')
axis tight;

% figure
subplot(224)
plot(covy_g1)
grid
title('Covarianza innovaciones y')
axis tight;

h1.Position=[0 0 1200 700];
h1.PaperUnits='points';
h1.PaperSize=[1200 700];
print('../Informe/Figuras/graf_ej3a','-dpdf','-bestfit');

% Grafico de medida, estimada, ruidosa
h2=figure;
subplot(221)
hold on
grid
plot(yk(:,1),yk(:,2),'color',myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x2(1,:),x2(2,:),'--b','LineWidth',2)
title('Estimación');
if(EsMatlab == 1)
    legend('Medida','Real','Estimada','location','SouthEast');
    xlabel('Posición x');
    ylabel('Posición y');
else
    legend(['Medida';'Real','Estimada'],'location','SouthEast');
    xlabel('Posición $x$ [\si{\m}]');
    ylabel('Posición $y$ [\si{\m}]');
end

% figure;
subplot(222)
plot(x'-x2');
xlim([0 150]);
title('Error de estimación');
xlabel('Tiempo');
legend('e_p_x','e_p_y','e_v_x','e_v_y','e_a_x','e_a_y');

covx_g2 = xcorr(g2(1,:)');
covy_g2 = xcorr(g2(2,:)');

% figure
subplot(223)
plot(covx_g2)
grid
title('Covarianza innovaciones x')
axis tight;

% figure
subplot(224)
plot(covy_g2)
grid
title('Covarianza innovaciones y')
axis tight;

h2.Position=[0 0 1200 700];
h2.PaperUnits='points';
h2.PaperSize=[1200 700];
print('../Informe/Figuras/graf_ej3b','-dpdf','-bestfit');

% Grafico de medida, estimada, ruidosa
h3=figure;
subplot(221)
hold on
grid
plot(yk(:,1),yk(:,2),'color',myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x3(1,:),x3(2,:),'--b','LineWidth',2)
title('Estimación');
if(EsMatlab == 1)
    legend('Medida','Real','Estimada','location','SouthEast');
    xlabel('Posición x');
    ylabel('Posición y');
else
    legend(['Medida';'Real','Estimada'],'location','SouthEast');
    xlabel('Posición $x$ [\si{\m}]');
    ylabel('Posición $y$ [\si{\m}]');
end

% figure;
subplot(222)
plot(x'-x3');
xlim([0 150]);
title('Error de estimación');
xlabel('Tiempo');
legend('e_p_x','e_p_y','e_v_x','e_v_y','e_a_x','e_a_y');

covx_g3 = xcorr(g3(1,:)');
covy_g3 = xcorr(g3(2,:)');

% figure
subplot(223)
plot(covx_g3)
grid
title('Covarianza innovaciones x')
axis tight;

% figure
subplot(224)
plot(covy_g3)
grid
title('Covarianza innovaciones y')
axis tight;

h3.Position=[0 0 1200 700];
h3.PaperUnits='points';
h3.PaperSize=[1200 700];
print('../Informe/Figuras/graf_ej3c','-dpdf','-bestfit');

% Grafico de medida, estimada, ruidosa
h4=figure;
subplot(221)
hold on
grid
plot(yk(:,1),yk(:,2),'color',myGreen)
plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
plot(x4(1,:),x4(2,:),'--b','LineWidth',2)
title('Estimación');
if(EsMatlab == 1)
    legend('Medida','Real','Estimada','location','SouthEast');
    xlabel('Posición x');
    ylabel('Posición y');
else
    legend(['Medida';'Real','Estimada'],'location','SouthEast');
    xlabel('Posición $x$ [\si{\m}]');
    ylabel('Posición $y$ [\si{\m}]');
end

% figure;
subplot(222)
plot(x'-x4');
xlim([0 150]);
title('Error de estimación');
xlabel('Tiempo');
legend('e_p_x','e_p_y','e_v_x','e_v_y','e_a_x','e_a_y');

covx_g4 = xcorr(g4(1,:)');
covy_g4 = xcorr(g4(2,:)');

% figure
subplot(223)
plot(covx_g4)
grid
title('Covarianza innovaciones x')
axis tight;

% figure
subplot(224)
plot(covy_g4)
grid
title('Covarianza innovaciones y')
axis tight;

h4.Position=[0 0 1200 700];
h4.PaperUnits='points';
h4.PaperSize=[1200 700];
print('../Informe/Figuras/graf_ej3d','-dpdf','-bestfit');

%%%%%%%%%% CÃ³mo grafico el error?!?! %%%%%%%%%%%%%%%%%%
