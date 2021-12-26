function t=inIR(x,y)
 % 2021-12-06 
##      if isscalar(x) && ~isscalar(y)
##         x = x .* ones(size(y));
##      endif
##      if ~isscalar(x) && isscalar(y)
##         y = y .* ones(size(x));
##      endif
##      if ~all(size(x.inf) == size(y.inf))
##        error("wrong dimensions");
##      endif
      t = inf(x) >= inf(y) & sup(x) <= sup(y);
endfunction 