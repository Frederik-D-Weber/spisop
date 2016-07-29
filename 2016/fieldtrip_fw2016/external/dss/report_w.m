function [report_data] = report_w(report_data, state)
% Reporting function to record w for each iteration
%

% Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
% Distributed by Laboratory of Computer and Information Science,
% Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.
% $Id: report_w.m,v 1.5 2005/04/20 10:19:24 kosti Exp $

report_data.w(:, state.iteration) = state.w;
