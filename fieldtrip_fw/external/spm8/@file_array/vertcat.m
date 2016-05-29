function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: vertcat.m 2781 2011-02-03 10:48:53Z roboos $


o = cat(1,varargin{:});
return;
