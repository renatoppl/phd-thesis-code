% Renato Paes Leme :: renatoppl@cs.cornell.edu
% this code implements an algorithm described in chapter 3 of my phd thesis:
%
% "Design and Analysis of Sponsored Search Mechanisms"
% Renato Paes Leme, PhD Thesis, Cornell University, 2012
%
% finds the cost of efficiency of the instance that has allocation
% pi and click-through-rates alpha. the cost of efficiency is the
% highest possible ratio between the highest revenue equilibrium and
% the highest ratio efficient equilibrium. Therefore, if there is
% an efficient revenue-optimal equilibrium, the cost of efficiency is 1.
% 
% assumptions:
% alpha(1) >= alpha(2) >= alpha(3) >= ...


function [beta, v, b] = cost_of_efficiency(alpha)
  beta = -1;
  n = length(alpha);
  for pi = perms([1:n])'
    for mask_number = 0:(2^n-1)
      mask = bitget(mask_number, [1:n]);
      [beta1, v1, b1] = cost_of_efficiency_sub(pi, mask, alpha);
      if (beta1 > beta)
        beta = beta1; v = v1; b = b1;
      endif
    endfor
  endfor
endfunction

% subroutine for computing the cost of efficiency
% assumptions:
% pi is an n-permutation, i.e., {pi(1), ..., pi(n)} = {1..n}  
% mask is an n-element {0,1}-vector

function [beta, v, b] = cost_of_efficiency_sub(pi, mask, alpha)
  % n = number of agents = number of slots
  n = length(alpha);

  % N = number of variables. We use var_v(i,n), var_b1(i,n) 
  % var_b2(i,n) and var_beta(n) to refer to variables 
  % v_i, b_i^1, b_i^2 and beta.
  N = 3*n+1;

  % naming of varibles
  var_v = inline("i", "i", "n");
  var_b1 = inline("i+n", "i", "n");
  var_b2 = inline("i+2*n", "i", "n");
  var_beta = inline("3*n+1", "n");

  % continuous variable types
  vartype = ""; for i=1:N vartype = strcat(vartype, "C"); endfor

  % bounds on the variables : between 0 and infty
  lb = zeros(N,1);
  ub = Inf * ones(N,1);

  % objective function : maximize the cost of efficiency stored as beta
  c = zeros(N,1); c(var_beta(n)) = 1;

  % direction of optimization = maximization (-1)
  sense = -1;

  % linear program constraints stored in A * x <= b
  % where A is an M * N matrix and b is an M vector.
  % const keeps track of the current constraint
  A = []; b = []; const = 1; ctype = "";

  % ordering of values constraints: v_i >= v_{i+1}, i=1..n-1
  for i=1:n-1
    A(const, var_v(i,n)) = 1;
    A(const, var_v(i+1,n)) = -1;
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor

  % ordering of bids constraints: b1_{pi(i)} >= b1_{pi(i+1)}, i=1..n-1
  for i=1:n-1
    A(const, var_b1(pi(i),n)) = 1;
    A(const, var_b1(pi(i+1),n)) = -1;
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor

  % non-overbidding constraints: b1_i <= v_i
  for i=1:n
    A(const, var_v(i,n)) = 1;
    A(const, var_b1(i,n)) = -1;
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor

  % Nash equilibrium constraints for k < j
  % alpha(j) * (v(pi(j)) - b1(pi(j+1))) >= alpha(k) * (v(pi(k)) - b1(pi(k)))
  for j=1:n for k=1:j-1
    A(const, var_v(pi(j),n)) = alpha(j) - alpha(k);
    if (j != n) A(const, var_b1(pi(j+1),n)) = - alpha(j); endif
    A(const, var_b1(pi(k),n)) = alpha(k); 
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor endfor

  % Nash equilibrium constraints for k < j
  % alpha(j) * (v(pi(j)) - b1(pi(j+1))) >= alpha(k) * (v(pi(k)) - b1(pi(k+1)))
  for j=1:n for k=j+1:n
    A(const, var_v(pi(j),n)) = alpha(j) - alpha(k);
    if (j != n) A(const, var_b1(pi(j+1),n)) = - alpha(j); endif
    if (k != n) A(const, var_b1(pi(k+1),n)) = alpha(k); endif 
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor endfor

  % we have that
  % b2(i) = min(v(i), ((alpha(i-1)-alpha(i)) * v(i-1) + alpha(i) *
  %                   b2(i+1))/alpha(i-1) )
  % given mask, we restrict our attention to (v,alpha) such that the left
  % is the minumum, i.e. b2(i) = v(i),  if mask(i) == 0 and b2(i) is defined by
  % the right if mask(i) == 1. So, for each alpha, we essentially partition the
  % space in 2^n parts (one for each mask) to takes care of the direction chosen
  % by the minimum.

  for i=2:n
    A(const, var_v(i,n)) = 1;
    A(const, var_v(i-1,n)) = - (alpha(i-1) - alpha(i)) / alpha(i-1);
    if (i<n) A(const, var_b2(i+1,n)) = - alpha(i) / alpha(i-1); endif
    b(const++) = 0;

    if (mask(i) == 0) 
      % v(i) >= ((alpha(i-1)-alpha(i))*v(i-1) + alpha(i) *b2(i+1))/alpha(i-1)
      % b2(i) == ((alpha(i-1)-alpha(i))*v(i-1) + alpha(i) *b2(i+1))/alpha(i-1)
      ctype = strcat(ctype, "L");
      A(const, var_v(i-1,n)) = - (alpha(i-1) - alpha(i)) / alpha(i-1);
      if (i<n) A(const, var_b2(i+1,n)) = - alpha(i) / alpha(i-1); endif
      A(const, var_b2(i,n)) = 1;
      b(const++) = 0;
      ctype = strcat(ctype, "S");
    elseif
      % v(i) <= ((alpha(i-1)-alpha(i))*v(i-1) + alpha(i) *b2(i+1))/alpha(i-1)
      % b2(i) == v(i)
      ctype = strcat(ctype, "U");
      A(const, var_v(i,n)) = -1;
      A(const, var_b2(i,n)) = 1;
      b(const++) = 0;
      ctype = strcat(ctype, "S");
    endif
  endfor


  % Normalization and defition of beta as the cost of efficiency
  for i=2:n  A(const, var_b2(i,n)) = alpha(i-1); endfor
  b(const++) = 1;
  ctype = strcat(ctype, "S");


  for i=2:n  A(const, var_b1(pi(i),n)) = alpha(i-1); endfor
  A(const, var_beta(n)) = -1;
  b(const++) = 0;
  ctype = strcat(ctype, "S");

  [xmin, fmin, status, extra] = ... 
    glpk(c, A, b, lb, ub, ctype, vartype, sense);

  % check for infeasibility (182,213 == primal-infeasible)
  if (status == 182 || status == 213)
    beta = -1;
    b = [];
    v = [];
  else
    b = [];
    beta = xmin(var_beta(n));
    for i=1:n
      v(i) = xmin(var_v(i,n));
      b(i) = xmin(var_b1(i,n));
    endfor
  endif
endfunction


