# !/bin/octave qf

function [bool_out]=set_graph(plot_function, LABELS, LEG, LOCATION, AXIS, bool_grid_minor)

	PLOTYY='plotyy';	%Cadenas comparativas
	PLOT='plot';
	lplot=length(PLOT);	%El script se encarga de comparar longitudes,  no las cadenas
	lplotyy=length(PLOTYY);
	lstyle=length(plot_function);


	if(bool_grid_minor)
		grid minor
	else			
		grid
	endif

	if(lstyle==lplotyy)
		xlabel(LABELS(1,:));
		ylabel(ax(1),LABELS(2,:));
		ylabel(ax(2),LABELS(3,:));
		
		axis(ax(1),AXIS(1,:));
		axis(ax(2),AXIS(2,:));
	else
		xlabel(LABELS(1,:),'FontSize', 14);
		ylabel(LABELS(2,:),'FontSize', 14);

		axis(AXIS);
	endif

		[Hleg, Hobj] = legend(gca,LEG,'location',LOCATION);
		legend('boxon');
		textobj = findobj(Hobj, 'type', 'text');
		set(textobj, 'Interpreter', 'LaTeX');
		set(textobj,'FontName','Arial','FontSize',14);

	bool_out=1;

endfunction
