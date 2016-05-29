function [varargout] = read_vista_mesh(varargin)

% READ_VISTA_MESH is implemented as mex file
%
% Use as
%   [nodes,elements,labels] = read_vista_mesh(filename);
% where
%   filename = the name of the Vista mesh (with extension .v)
%
% $Id: read_vista_mesh.m 7197 2012-12-15 10:49:56Z roboos $

error('The mex file %s is missing', [mfilename '.' mexext]);
