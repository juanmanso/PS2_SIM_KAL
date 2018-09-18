
function [bool]=validate_arguments(arg_vector,argc)
	bool=0;	
	
	if(argc!=4)
		printf("\n\tERROR: Cantidad de argumentos inválida.\n\tInvocar de la siguiente manera:\n\t\t$ ./tp3 -p <bool> -results <bool>\n");
		return;
	elseif(arg_vector{1}!="-p")
		printf("\n\tERROR: Argumento 1 inválido.\n\tInvocar de la siguiente manera:\n\t\t$ ./tp3 -p <bool> -results <bool>\n");
		return;
	elseif(str2num(arg_vector{2})!=1 && str2num(arg_vector{2}))
		printf("\n\tERROR: Argumento 2 booleano inválido.\n\tInvocar de la siguiente manera:\n\t\t$ ./tp3 -p <bool> -results <bool>\n");
		return;
	elseif(arg_vector{3}!="-results")
		printf("\n\tERROR: Argumento 3 inválido.\n\tInvocar de la siguiente manera:\n\t\t$ ./tp3 -p <bool> -results <bool>\n");
		return;
	elseif(str2num(arg_vector{4})==1 && !str2num(arg_vector{4}))
		printf("\n\tERROR: Argumento 4 booleano inválido.\n\tInvocar de la siguiente manera:\n\t\t$ ./tp3 -p <bool> -results <bool>\n");
		return;	
	endif

	bool=1;

endfunction

