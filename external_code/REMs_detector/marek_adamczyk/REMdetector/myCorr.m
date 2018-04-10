function rho = myCorr(x,y)
%CORR     Correlation coefficient.
%
%         rho = corr(x,y)
%
%         The correlation coeffient is computed between all
%	  columns in  x  and all coluns in  y. If only one
%	  argument  x  is given then the correlation matrix
%	  for  x  is returned instead.
%
%	  See also CVAR, SPEARMAN.

%       GPL Copyright (c) Anders Holtsberg, 1998

if nargin < 2
   c = myCvar(x);
   s = sqrt(diag(c));
   rho = c ./ (s*s');
   n = size(x,2);
   rho(1:n+1:n^2) = ones(n,1);
else
   cx = myCvar(x);
   sx = sqrt(diag(cx));
   cy = myCvar(y);
   sy = sqrt(diag(cy));
   c = myCvar(x,y);
   rho = c ./ (sx*sy');
end
