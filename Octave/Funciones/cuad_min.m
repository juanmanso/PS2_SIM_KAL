


function [] = cuad_min(M,w0,a,W)

	if(mod(M,2)==0) % Si M es par
		Fw = @F_tipo_I;
		k = M/2;
	
	else
		Fw = @F_tipo_II;
		k = (M-1)/2;
	end

	w = (0:0.01:1)*pi;
	m = 0:k;
	b = m*0; %Inicializo b
	c = b;	 %Inicializo c

	for i=1:length(W)
		b += 2*W(i)*(feval(@Fw,w)').^2.*cos(w'.*m).*(w>w0((2*i)-1) && w<w0(2*i));
		c += 2*W(i)*a(i)*(feval(@Fw,w)').*cos(w'.*m).*(w>w0((2*i)-1) && w<w0(2*i));
	end

	g = c^(-1)*b;

end
