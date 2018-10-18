config_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EJ KALMAN - Perdida de datos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datos_str = load('datos.mat');

Acel = datos_str.Acel;
Tiempo = datos_str.tiempo;
Pos = datos_str.Pos;
Vel = datos_str.Vel;

dim = 2;			% Se considera solo x e y
tipos_variables = 3;		% Posicion, Velocidad, Aceleracion
cant_mediciones = length(Pos);
cant_estados = tipos_variables * dim;

%% p_bernoulli = 0.90;		% Prob de 1
%% p_perdido = 1 - p_bernoulli;	% Prob de 0
%p_perdido_p = rand(1,dim)*val_max_perdido + (1-val_max_perdido);
%p_perdido_v = rand(1,dim)/2 + 0.5; 
%p_perdido_a = rand(1,dim)/2 + 0.5;  

% Variables de configuracion de la perdida de datos
p_perdido = [0.1 0.9];				% Se pierden 3 de 10 datos
for i=1:length(p_perdido)
    p_perdidoaux(i) = 1-p_perdido(i);
    p_perdido_p(i,:) = [1 1].*p_perdidoaux(i);
    p_perdido_v(i,:) = [1 1].*p_perdidoaux(i);
    p_perdido_a(i,:) = [1 1].*p_perdidoaux(i);
    p_bernoulli(i,:) = [p_perdido_p(i,:), p_perdido_v(i,:), p_perdido_a(i,:)];
end
% Hay mediciones de:
bool_p = 1;
bool_v = 0;
bool_a = 0;

% Seleccion de impresion de imagenes
bool_print = 1;

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
[rowC,colC]=size(C);

M_eta = [randn(dim,cant_mediciones)*sigma_etap*bool_p; 
	randn(dim,cant_mediciones)*sigma_etav*bool_v;
       	randn(dim,cant_mediciones)*sigma_etaa*bool_a];

Mediciones = [Pos(:,1:dim) Vel(:,1:dim) Acel(:,1:dim)];
%Mediciones = Mediciones .* binornd(1,1-p_bernoulli,[cant_estados,cant_mediciones])';

yk = C * Mediciones' + (C*M_eta);
yk = yk'; % Asi tiene la forma de Pos
%yk = yk' .* binornd(1,p_bernoulli,[rowC,cant_mediciones])'; % Asi tiene la forma de Pos

R = diag([ones(1,dim*bool_p)*sigma_etap^2 ones(1,dim*bool_v)*sigma_etav^2 ones(1,dim*bool_a)*sigma_etaa^2]);


for l=1:length(p_perdido)
    
    %%% ALGORITMO %%%%
    x = x0;
    P = P0_0;
    xk1_k1 = x;
    Pk1_k1 = P;
    g = yk(1,:)';

    for i=1:cant_mediciones-1
        % Perdida de datos
        perder_dato = binornd(1,p_bernoulli(l,:),[1,cant_estados]);
        
        C_aux = C;
        if EsMatlab
            C=bsxfun(@times,C,perder_dato);
        else % Octave
            C = C .* perder_dato;
        end
        
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
        C = C_aux;
        
        
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
    
    suptitle(sprintf('p_{perdida}=%.2f',p_perdido(l)));
    
    h.Position=[0 0 1200 700];
    h.PaperUnits='points';
    h.PaperSize=[1200 700];
    if bool_print
        print(sprintf(['../Informe/Figuras/graf_ej7_' num2str(l)]),'-dpdf','-bestfit');
    end
    
    % Grafico de correlacion de innovaciones (debe ser ruido blanco)
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
    
    suptitle(sprintf('p_{perdida}=%.2f',p_perdido(l)));
    
    h2.Position=[0 0 1200 700];
    h2.PaperUnits='points';
    h2.PaperSize=[1200 700];
    if bool_print
        print(sprintf(['../Informe/Figuras/covinn_ej7_' num2str(l)]),'-dpdf','-bestfit');
    end
end

% % Grafico de medida, estimada, ruidosa
% figure
% hold on
% grid
% plot(x(1,:),x(2,:),'LineWidth',3)
% plot(Pos(:,1),Pos(:,2),'r','LineWidth',2)
% plot(yk(:,1),yk(:,2),'color',myGreen)
% title('Estimacion');
% if(EsMatlab == 1)
%     legend('Estimada','Medida','Ruidosa');
% else
%     legend(['Estimada';'Medida';'Ruidosa']);
% end
% xlabel = 'Tiempo [s]';
% ylabel = 'Posicion [m]';
% 
% % Grafico del estado posicion en funcion del tiempo
% figure
% hold on
% grid
% plot(x(1,:),'LineWidth',2)
% plot(x(2,:),'color',myGreen,'LineWidth',2)
% title('Estados de posicion');
% 
% 
% % Grafico del estado velocidad en funcion del tiempo
% figure
% hold on
% grid
% plot(x(3,:),'LineWidth',2)
% plot(x(4,:),'color',myGreen,'LineWidth',2)
% title('Estados de velocidad');
% 
% 
% % Grafico del estado aceleracion en funcion del tiempo
% figure
% hold on
% grid
% plot(x(5,:),'LineWidth',2)
% plot(x(6,:),'color',myGreen,'LineWidth',2)
% title('Estados de aceleración');
% 
% 
% % Grafico de correlacion de innovaciones (debe ser ruido blanco)
% covx_g = xcorr(g(1,:)');
% covy_g = xcorr(g(2,:)');
% 
% figure
% plot(covx_g)
% grid
% title('Covarianza innovaciones x')
% 
% figure
% plot(covy_g)
% grid
% title('Covarianza innovaciones y')




% Observabilidad
Obs = obsv(Ad,C);
rango_obs = rank(Obs);
estados_no_observables = cant_estados - rango_obs


