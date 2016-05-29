function o = horzcat(varargin)
% Horizontal concatenation of file_array objects
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

%
% $Id: horzcat.m 2781 2011-02-03 10:48:53Z roboos $

o    = cat(2,varargin{:});
return;

