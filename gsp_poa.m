% Renato Paes Leme :: renatoppl@cs.cornell.edu
% this code implements an algorithm described in chapter 2 of my phd thesis:
%
% "Design and Analysis of Sponsored Search Mechanisms"
% Renato Paes Leme, PhD Thesis, Cornell University, 2012
%
% finds the price of anarchy of the instance that has allocation
% pi and click-through-rates alpha.
% 
% assumptions:
% length(pi) == length(alpha)
% alpha(1) >= alpha(2) >= alpha(3) >= ...
% pi is a permutation, i.e., {pi(1), ..., pi(n)} = {1, ... , n}

function [beta, v, b] = gsp_poa(pi, alpha)
  % n = number of agents = number of slots
  n = length(alpha);

  % N = number of variables. We use var_v(i,n), var_b(i,n) 
  % and var_beta(n) to refer to variables v_i, b_i and beta.
  % see functions below that define the variable indexing
  N = 2*n+1;

  % naming of varibles
  var_v = inline("i", "i", "n");
  var_b = inline("i+n", "i", "n");
  var_beta = inline("2*n+1", "n");

  % continuous variable types
  vartype = ""; for i=1:N vartype = strcat(vartype, "C"); endfor

  % bounds on the variables : between 0 and infty
  lb = zeros(N,1);
  ub = Inf * ones(N,1);

  % objective function : maximize price of anarchy
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

  % ordering of bids constraints: b_{pi(i)} >= b_{pi(i+1)}, i=1..n-1
  for i=1:n-1
    A(const, var_b(pi(i),n)) = 1;
    A(const, var_b(pi(i+1),n)) = -1;
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor

  % non-overbidding constraints: b_i <= v_i
  for i=1:n
    A(const, var_v(i,n)) = 1;
    A(const, var_b(i,n)) = -1;
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor

  % Nash equilibrium constraints for k < j
  % alpha(j) * (v(pi(j)) - b(pi(j+1))) >= alpha(k) * (v(pi(k)) - b(pi(k)))
  for j=1:n for k=1:j-1
    A(const, var_v(pi(j),n)) = alpha(j) - alpha(k);
    if (j != n) A(const, var_b(pi(j+1),n)) = - alpha(j); endif
    A(const, var_b(pi(k),n)) = alpha(k); 
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor endfor

  % Nash equilibrium constraints for k < j
  % alpha(j) * (v(pi(j)) - b(pi(j+1))) >= alpha(k) * (v(pi(k)) - b(pi(k+1)))
  for j=1:n for k=j+1:n
    A(const, var_v(pi(j),n)) = alpha(j) - alpha(k);
    if (j != n) A(const, var_b(pi(j+1),n)) = - alpha(j); endif
    if (k != n) A(const, var_b(pi(k+1),n)) = alpha(k); endif 
    b(const++) = 0;
    ctype = strcat(ctype, "L");
  endfor endfor

  % Normalization and Price of Anarchy constraints
  for i=1:n  A(const, var_v(i,n)) = alpha(i); endfor
  A(const, var_beta(n)) = -1;
  b(const++) = 0;
  ctype = strcat(ctype, "S");


  for i=1:n  A(const, var_v(pi(i),n)) = alpha(i); endfor
  b(const++) = 1;
  ctype = strcat(ctype, "S");

  [xmin, fmin, status, extra] = ... 
    glpk(c, A, b, lb, ub, ctype, vartype, sense);

  b = [];
  beta = xmin(var_beta(n));
  for i=1:n
    v(i) = xmin(var_v(i,n));
    b(i) = xmin(var_b(i,n));
  endfor
endfunction


