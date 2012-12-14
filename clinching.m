% Renato Paes Leme :: renatoppl@cs.cornell.edu
% this code implements an algorithm described in chapter 5 of my phd thesis:
%
% "Design and Analysis of Sponsored Search Mechanisms"
% Renato Paes Leme, PhD Thesis, Cornell University, 2012
%
% computes the allocation and payments in the Adaptive Clinching Auction
% given the values of the players v (which is a vector of n components),
% the budgets B of the players and the supply.
%
% assumptions:
% v and B are vectors of the same size (call it n)
% s is a positive scalar
%
% the output consists of two vectors x and pi of size n such that
% \sum_{i=1}^n x(i) = s(i) and pi(i) <= B(i) for all i


function [x,pi] = clinching(v,B,s, tol=1E-6)

  % store the vector of initial budgets
  B0 = B;

  % number of agents
  n = length(v); 
  
  % allocation : x(i) = amount of goods player i gets
  x = zeros(1,n);

  % mantains the set of active players, i.e., A = {i; d(i) > 0} 
  % at the current price
  A = [1:n];

  % maintains the clinching set. C is the set of players clinching
  C = [];

  % initial price (a small enough price that no clinching happens
  % before it).
  p = .5 * min( min(v), (sum(B) - max(B)) / s);

  not_finished = 1;

  while (not_finished)
    % we consider two types of events: either the price clock reaches the
    % value of a certain player and it drops out of the auction (event 1)
    % or, one extra player starts clinching (event 2).

    % event 1 : stores in E1 the value of the active player with minimum
    % value and in iE1 the index of such player
    [E1, iE1] = min(v(A));
    iE1 = A(iE1);

    % event 2
    E2 = compute_event_2(s,p,B,A,C);
  
    p
    E1
    E2
    str= "--------------------------"

    if (E1 <= E2)
      % calculates allocation and payment updates since last event (price
      % p), i.e., amount clinched by players until E1- (just before 
      % price reaches E1)
      c = length(C);
      delta = zeros(n,1);
      for i = C
        delta(i) = (s/c)*(1-((p/E1)**c));
        if c==1
          B(i) -= s*p*(log(E1) - log(p));
        else
          B(i) -= (s*p/(c-1))*(1-((p/E1)**(c-1)));
        endif
      endfor
      for i=C
        x(i) += delta(i);
        s -= delta(i);
      endfor

      r=1
      delta
      
      % now we see how allocation and payments get updated as the price
      % reaches E1, i.e., clinching that happens between E1- and E1.
      p = E1;

      % removes player iE1 from the active set
      A = setdiff(A,iE1);
      if ismember(iE1, C)
        % the auction finishes once one player in the clinching set drops
        not_finished = 0;
      endif

      % for the remaining players, check how much they clinch at
      % price p = E1. If they clinched a positive amount and wasn't in the
      % clinching set before, we add them to the clinching set.
      delta = zeros(n,1);
      for i=A
        if ((sum(B(A)) - B(i))/p < s)
          if !ismember(i,C)
            C = union(C,[i]);
          endif
          delta(i) = s - (sum(B(A)) - B(i))/p;
        endif
      endfor
      for i=C
        x(i) += delta(i);
        B(i) -= p*delta(i);
        s -= delta(i);
      endfor
      if length(A) == 1
        not_finished = 0;
      endif
    else
      % in the case event 2 happens, i.e., a new player enters
      % the clinching set, first we calculate the player with largest
      % budget who is not in the clinching set.
      Bnc = max(B(setdiff(A,C)));
      delta = zeros(n,1);
      c = length(C);
      for i=C
        delta(i) = (s/c)*(1-((p/E2)**c));
      endfor
      for i=C
        x(i) += delta(i);
        s -= delta(i);
        B(i) = Bnc;
      endfor
      for i = setdiff(A,C)
        if abs(B(i) - Bnc) < tol
          C  = union(C,[i]);
        endif
      endfor

      r=2
      p
      E1
      E2
      delta

      p = E2;

    endif
  endwhile

  % the payment is the difference between the original budget and the
  % remaining budget
  pi = B0 - B;
endfunction


% event 2 corresponds to the next time a player will enter the clinching set
function E2 = compute_event_2(s,p,B,A,C)
  if length(C) == length(A)
    % if all active players are in the clinching set, no additional player will
    % enter the clinching set, so the event happens for price infinity
    E2 = Inf;
  elseif length(C) == 0
    % if there is no player clinching yet, the player with the highest budget
    % enters the clinching set in the moment the aggregate demand of all the
    % other players equals the supply, ( sum(B(A)) - max(B(A)) ) / p = s
    E2 = ( sum(B(A))-max(B(A)) ) / s;
  else
    % we define Bc to be the "clinching budget", which is the value of the
    % budget in the clinching set, i.e., max(B(A)) = B(i) for any i in C.
    % and BnC the largest "non-clinching budget" which is max(B(A \setminus C))
    % then we calculate E2 as being the price for which Bc will come down
    % to the value of Bnc
    Bc = max(B(C));
    Bnc = max(B(setdiff(A,C)));
    c = length(C)
    if (c == 1)
      E2 = exp((Bc-Bnc + s*p* log(p))/(s*p));
    else
      E2 = ( (1/(p**(c-1))) + (Bnc - Bc)*(c-1) / (s * p**c) )**(-1/(c-1));
    endif
  endif
endfunction
