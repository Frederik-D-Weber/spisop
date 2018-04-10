function c = myCvar(x,y)
%CVAR     Covariance.
%
%         c = cvar(x,y)
%
%         The covariance is computed between all columns in  x  and 
%	  all columns in  y. If only one argument  x  is given then 
%	  the covariance matrix for  x  is returned instead, i e
%	  the result is the same as for the call  cvar(x,x).
%
%	  See also CORR, SPEARMAN.

%       GPL Copyright (c) Anders Holtsberg, 1998, 1999

if nargin < 2
   if size(x,1) == 1
      x = x';
   end
   n = size(x,1);
   m = mean(x);
   x = x-m(ones(n,1),:);
   c = (x'*x)/(n-1);
else
   if size(x,1) == 1
      x = x';
   end
   if size(y,1) == 1
      y = y';
   end
   n = size(x,1);
   mx = mean(x);
   my = mean(y);
   x = x-mx(ones(n,1),:);
   y = y-my(ones(n,1),:);
   c = (x'*y)/(n-1);   
end
