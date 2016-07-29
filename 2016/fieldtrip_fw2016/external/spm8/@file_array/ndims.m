function out = ndims(fa)
% Number of dimensions
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: ndims.m 2781 2011-02-03 10:48:53Z roboos $


out = size(fa);
out = length(out);

