function tf = isfield(this,field)
% Isfield method for GIfTI objects
% FORMAT tf = isfield(this,field)
% this   -  GIfTI object
% field  -  string of cell array
% tf     -  logical array
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: isfield.m 3261 2011-03-31 15:06:56Z roboos $

tf = ismember(field, fieldnames(this));