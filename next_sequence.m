% Renato Paes Leme :: renatoppl@cs.cornell.edu
% this code implements is an auxiliary function used by the other functions
% in this folder
%
% provides a "next"-function to generate all increasing sequences of 1..k
% of n elements. 
%
% usage:
%s = ones(1,n);
%do
%  ... do something with s ...
%  s = next_sequence(s, k);
%until (length(s) == 0) 


function s = next_sequence(s, k)
  n = length(s);
  i = n;
  while (i > 0 && s(i) == k)
    --i;
  endwhile
  if (i == 0)
    s = [];
    return;
  endif
  s(i) += 1;
  for j=i+1:n
    s(j) = s(i);
  endfor
endfunction


