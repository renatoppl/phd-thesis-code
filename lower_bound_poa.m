% Renato Paes Leme :: renatoppl@cs.cornell.edu
% this code implements an algorithm described in chapter 2 of my phd thesis:
%
% "Design and Analysis of Sponsored Search Mechanisms"
% Renato Paes Leme, PhD Thesis, Cornell University, 2012
%
% this code uses calculates the worse Price of Anarchy (PoA) for GSP among all
% instances with n players and n slots, where valuations and budgets are between
% zero and one and multiples of epsilon. It does so by enumerating over all
% instances and calling the gsp_poa function for each such instance.

function [poa, alpha, v, b, pi] = lower_bound_poa(n, epsilon)
  poa = -1;
  alpha_vals = [1:-epsilon:0];
  s = ones(1,n-1);
  do 
    alpha1 = [1, alpha_vals(s)];
    s = next_sequence(s, length(alpha_vals));
    for pi1 = perms([1:n])'
      [beta, v1, b1] = gsp_poa(pi1, alpha1);
      if (beta > poa) 
        poa = beta; alpha = alpha1;
        v = v1; b = b1; pi = pi1';
      endif
    endfor
  until (length(s) == 0)  
endfunction
