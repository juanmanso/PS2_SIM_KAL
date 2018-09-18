
% Función auxiliar para imprimir "más rápido" las figuras

function [] = imprimir_figura(nombre)

	format_out_name=".tex";			% Configuro el formato de impresión
	format_print_config="-depslatex";	% Idem
	print_path="./graf_";	% Ruta de salida

	print([print_path,nombre,format_out_name],format_print_config, "-color");

endfunction
