function l = length(x)
% Overloaded length function for file_array objects
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: length.m 2781 2011-02-03 10:48:53Z roboos $


l = max(size(x));

