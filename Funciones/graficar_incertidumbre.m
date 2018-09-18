

function graficar_incertidumbre(P,t,a,s)
	PA = P(1,:); PB = P(2,:); PC = P(3,:); PD = P(4,:);
	sbx =s(1); sby=s(2); ssx=s(3); ssy=s(4); srx=s(5); sry=s(6);
	tf = t(end);

	dpx = 3*(tf^2/2 *(sbx+max(abs(a(:,1)))*ssx) + srx) 
	dpy = 3*(tf^2/2 *(sby+max(abs(a(:,2)))*ssy) + sry) 
	% Da muy chico me parece... no encuentro el error

	h = linspace(0,2*pi,10);

	hold on
	plot(dpx*cos(h)+P(1,1), dpy*sin(h)+P(1,2),'r.-')
	plot(dpx*cos(h)+P(2,1), dpy*sin(h)+P(2,2),'r.-')
	plot(dpx*cos(h)+P(3,1), dpy*sin(h)+P(3,2),'r.-')
	plot(dpx*cos(h)+P(4,1), dpy*sin(h)+P(4,2),'r.-')
endfunction
