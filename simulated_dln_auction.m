function [x,p] = simulated_dln_auction(v,B,T)
  n = length(v);
  x = zeros(n,1); 
  p = zeros(n,1);
	r = 1.0; %remnant supply
	eps = 1. / double(T);
	pc = zeros(n,1);
	for j=1:n
		pc(j) = eps;
  endfor
	ii = 1;
	while (pc(n) < max(v)+eps)
		pc(ii) += eps;
		ii = 1 + mod(ii,n);
		%pc
		%B
		%transpose(x)
		%transpose(p)
		%calculate demands
		d = []; D = 0;
		for i=1:n
			if pc(i) <= v(i)
				d(i) = double(B(i)) / pc(i);
			else
				d(i) = 0;
			endif
			D += d(i);
		endfor
		%d
		%r
		%calculate clinched amounts
		clinch = [];
		for i=1:n
			clinch(i) = max(0,r - (D - d(i)));
		endfor
		for i=1:n
			x(i) += clinch(i);
			r -= clinch(i);
			p(i) += pc(i) * clinch(i);
			B(i) -= pc(i) * clinch(i);
			B(i) = max(B(i), 0);
		endfor
		%clinch
		%disp("----------------------------\n")
	endwhile
	x = transpose(x);
	p = transpose(p);
endfunction
