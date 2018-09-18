# Escrito por Juan Manso, Facultad de Ingeniería UBA.

# La siguiente función se encarga de graficar en estilo plotyy o plot
# Los vectores a graficar se pasan como columnas en MATRIX
# Se supone que el vector "tiempo" es la primera columna
# Las columnas subsiguientes son valores de funciones para dichos "t".
# Si se quiere graficar con plotyy, entonces la matriz tiene que tener 3 columnas
# Si se quiere graficar con plot, no importa las columnas, las grafica todas con el mismo tipo de linea
# Si se quiere graficar con errorbar, con los errores constantes para todos los puntos, se ingresa con "errorbart" y los errores son las columnas de AUX.

% MATRIX tiene en las columnas las variables y en las filas los valores de cada Var para un "t".
% style es una cadena de caracteres con el nombre de la función a utilizar
% OUTPUT_FILE_NAME es una cadena de caracteres que tiene el nombre y ubicación del archivo
% AUX es una variable que lleva el contenido extra necesario para manejar ciertas funciones como por ejemplo errbar.
% LABELS es un arreglo de cadenas de caracteres que contiene las etiquetas para los ejes
% LEG es un arreglo de cadenas de caracteres que contiene el texto que se expone en las leyendas
% fig_number es un número. Si es 0, genera una nueva figura. Si es !=0, sobreescribe/ genera la figure(number)
% line_style es una cadena de caracteres con el tipo de linea en el gráfico [p.e. "o", "-","^", etc.]
% bool_grid_minor es un booleano. Si es 1, grafica con grid minor. Si es 0, grafica con grid común.

function [bool_out]=graph_same_length(MATRIX, AUX, style, OUTPUT_FILE_NAME, LABELS, LEG, LOCATION, AXIS, fig_number,line_style, bool_grid_minor)

	SEMILOG="semilog";
	PLOTYY="plotyy";	%Cadenas comparativas
	PLOT="plot";
	ERRORBAR="errorbar";
	lsemilog=length(SEMILOG);
	lplot=length(PLOT);	%El script se encarga de comparar longitudes,  no las cadenas
	lplotyy=length(PLOTYY);
	lerrorbar=length(ERRORBAR);
	lstyle=length(style);
	c=columns(MATRIX);	%Calculo la cantidad de columnas para ver como fue ingresada MATRIX


	% Configuro al programa para graficar como GNUPLOT (en vez de FLTK definido por defecto)
	graphics_toolkit("gnuplot");
	hold on;

	if(length(MATRIX)>c)	%Transpongo MATRIX de ser necesario
		M=MATRIX;
	else
		M=MATRIX';
		c=columns(M);
	endif


	MAX_COLUMS_yy=3;	%Macro para validar si se puede usar plotyy

	bool_out=0;
	X=1;
	Y=2;
	Z=3;


	if(fig_number) 		
		bool_new_figure=0;
	else
		bool_new_figure=1;
	endif

	if(lstyle!=lsemilog && lstyle!=lplot && lstyle!=lplotyy && lstyle!=lerrorbar)	%Validacion del estilo
		disp('Estilo insertado inválido')
		bool_out=0;
		return;
	endif

	if(isprint(OUTPUT_FILE_NAME))		%Si se ingresa un caracter no imprimible, no imprime el gráfico
		bool_print=0;
	else
		bool_print=1;
	endif
	
	
	
	if(!bool_new_figure)	%Si el numero de figura es nulo, que genere una nueva. Si no sobreescribe/crea una con ese numero
		figure(fig_number)
	endif


			

		if(lstyle==lplotyy)
			if(c!=MAX_COLUMS_yy)
				disp('Cantidad de vectores a graficar con plotyy inválido')
				bool_out=0;
				return;
			else
			ax=plotyy_style(M,X,Y,Z);
			endif
		elseif(lstyle==lerrorbar)
		
			if(!errorbar_style(M,X,AUX,line_style))
				bool_out=0;
				return;
			endif
		else
			plot_style(M,X,line_style);
		endif
		

	if(bool_grid_minor)
		grid minor
	else			
		grid
	endif
		

	if(lstyle==lplotyy)
		xlabel(LABELS(X,:));
		ylabel(ax(1),LABELS(Y,:));
		ylabel(ax(2),LABELS(Z,:));
		
		axis(ax(1),AXIS(1,:));
		axis(ax(2),AXIS(2,:));
	else
		xlabel(LABELS(X,:));
		ylabel(LABELS(Y,:));
		
		axis(AXIS);
	endif

		Hleg = legend(gca,LEG,'location',LOCATION);
		legend('boxon');
		set(Hleg,'FontName','Arial','FontSize',8);
	
	
	if(bool_print)
	n=1;
	n1=1;
	n2=1;
	check=0;

		while(n<=length(OUTPUT_FILE_NAME))	% Busca las carpetas necesarias para imprimir de ser necesario.
			if(OUTPUT_FILE_NAME(n)=='/')
				check++;
			elseif(OUTPUT_FILE_NAME=='.')
				n2=n;
			endif
			if(check==2)
				n1=n;
			endif	
			n++;
		endwhile
		if(n1!=1)
			mkdir ('./',OUTPUT_FILE_NAME(1:n1));
		endif
		if(n2!=1)
			if(OUTPUT_FILE_NAME(end-n2:n2) == "tex") % Si se quiere imprimir para LaTeX
				PRINT_COLOR="-color";    
				FORMAT_PRINT_CONFIG="-depslatex"; 		%Idem
				PRINT_PATH="./";                		% Ruta de salida

				if(n1!=1)
					print([OUTPUT_FILE_NAME],FORMAT_PRINT_CONFIG, PRINT_COLOR); 
				else
					print([PRINT_PATH,OUTPUT_FILE_NAME],FORMAT_PRINT_CONFIG, PRINT_COLOR); 
				endif
			else
				print(OUTPUT_FILE_NAME,['-d' OUTPUT_FILE_NAME(end-n2:n2)]);
			endif
		endif
	endif

	bool_out=1;
endfunction



function [ax]=plotyy_style(M,X,Y,Z)
	ax=plotyy(M(:,X),M(:,Y),M(:,X),M(:,Z));
endfunction

function []=plot_style(M,X,line_style)
	k=X;
	while(k<=rows(line_style) && (k+1<=columns(M)))
		plot(M(:,X),M(:,k+1),line_style(k,:));
		k++;
	endwhile
endfunction

function [bool]=errorbar_style(M,X,ERR,line_style);
	bool=1;
	hold on;
	j=1; % Cantidad de graficos
	k=1; % Auxiliar para recorrer las matrices
	while(k<=rows(line_style) && (k+1<=columns(M)))
		if(columns(ERR)>=(columns(M)-1)) % Si tiene error en y 
			if(columns(ERR)==2*(columns(M)-1)) % Si tiene error en x
				errorbar(M(:,X),M(:,k+1),ERR(:,j),ERR(:,j+1),strcat('~>', line_style(k,:)));

				j++;
			else
				errorbar(M(:,X),M(:,k+1),ERR(:,j),strcat('~',line_style(k,:)));
			endif
		else
			disp('Matriz de errores inválida para graficar con errorbar');
			bool=0;
			return;
		endif

		j++;
		k++;
	endwhile
endfunction
