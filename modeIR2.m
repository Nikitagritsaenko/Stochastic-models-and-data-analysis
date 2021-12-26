% Created: 2021-11-23
function [mode, mu_array, max_mu, mode_ind, c_array, C, multi]= modeIR2(X)
  % INI
  mu_array=[]; mode_ind=[]; c_array=[]; C=[];  multi=0;
  %
  % X=PgammaStdNow
  % X=DataStdCover
  % X=DataStd
  I=intersect(X(:));
  if not(isempty(I)) 
    mode=I, max_mu=length(X); multi=1;
    C=[inf(X)', sup(X)'];
    C=sort(C);
    % 2021-11-30 remove duplicate
    C=unique(C);
    for ii=1:length(C)-1
      c(ii)=infsupdec(C(ii),C(ii+1));
    end
    c_array=c;
    %
    %for ii=1:length(C)-1
    %mu(ii)=0;
    %  for jj=1:length(X)
    %      if not(isempty(intersect(c(ii), X(jj)))) 
    %        mu(ii)=mu(ii)+1;
    %      end
    %  end
    %end
    %mu_array=mu;
    %mode_ind=[];
    %for ii=1:length(c)
    %  if mu_array(ii)== max_mu
    %    mode_ind=[mode_ind, ii];
    %  endif
    %end
    %
  else   
    C = [ inf(X)', sup(X)' ];
    C = sort(C);
    % 2021-11-30 remove duplicate
    C = unique(C);
    for ii=1:length(C)-1
      c(ii)=infsupdec(C(ii),C(ii+1));
    end
    c_array=c;
    clear mu
    length(C)
    for ii=1:length(C)-1
      mu(ii)=0;
      ii
      for jj=1:length(X)
        %      if not(isempty(intersect(c(ii), X(jj)))) 
        if inIR(c(ii), X(jj)) > 0 
          mu(ii)=mu(ii)+1;
        end
      end
    end
    mu_array=mu;
    [max_mu,max_mu_ind]=max(mu_array);
    mode=[];
    mode_ind=[];
    multi=0;
    for ii=1:length(c)
      if mu_array(ii)== max_mu
        multi=multi+1
        mode=[mode, c(ii)];
        mode_ind=[mode_ind, ii];
      endif
    end
  #endif
endfunction
