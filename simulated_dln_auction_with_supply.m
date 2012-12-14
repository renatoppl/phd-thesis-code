function [x,p] = simulated_dln_auction_with_supply(v,B,T,s)
  [x,p] = simulated_dln_auction(v*s,B,T);
  x=x*s;
endfunction
