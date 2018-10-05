

function [gk] = innovaciones(yk, Ck, xk_k1)

	gk = yk' - Ck*xk_k1;

end

