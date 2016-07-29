% FT_PREAMBLE_PROVENANCE is a helper script that records the time and memory at the
% start of the function. This is to be used together with FT_POSTAMBLE_PROVENANCE which
% will record and store the time and memory at the end of the function. This is
% stored in the output configuration together with information about the enbvironment,
% such as the name of the user and computer, the matlab and fieldtrip version, etc.
%
% Another aspects of provenance relates to uniquely identifying the input and the
% output data. The code that deals with tracking the information about the input data
% structures is found in ft_preamble_loadvar. The code that deals with tracking the
% information about the output data structures is found in ft_preamble_history.

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: ft_preamble_provenance.m 8070 2013-04-24 13:38:40Z roevdmei $

% Record the start time and memory. These are used by ft_postamble_callinfo, which
% stores them in the output cfg.callinfo.  In the mean time, they are stored in the
% function workspace, which is why they should have cryptical names to prevent any
% variable name clashes.

% the name of the variables are passed in the preamble field
global ft_default

if (isfield(cfg, 'trackcallinfo') && ~istrue(cfg.trackcallinfo)) 
  % do not track the call information
  return
end

% add the user-specified cfg (before any defaults handling etc.) to the
% callinfo
cfg.callinfo.usercfg = cfg;

% compute the MD5 hash of each of the input arguments
% temporarily remove the cfg field for getting the hash (creating a duplicate of the data, but still has the same mem ref, so no extra mem needed)
if isequal(ft_default.preamble, {'varargin'})
  tmpargin = varargin;
else
  tmpargin = cellfun(@eval, ft_default.preamble, 'UniformOutput', false);
end
cfg.callinfo.inputhash = cell(1,numel(tmpargin));
for iargin = 1:numel(tmpargin)
  tmparg = tmpargin{iargin}; % can't get number of bytes with whos unless taken out of it's cell
  if isfield(tmparg,'cfg')
    tmparg = rmfield(tmparg,'cfg');
  else
  end
  % only calculate md5 when below 2^31 bytes (CalcMD5 can't handle larger input)
  bytenum = whos('tmparg');
  bytenum = bytenum.bytes;
  if bytenum<2^31
    cfg.callinfo.inputhash{iargin} = CalcMD5(mxSerialize(tmparg));
  end
end
clear tmpargin tmparg; % remove the extra references


stack = dbstack('-completenames');
% stack(1) is this script
% stack(2) is the calling ft_postamble function
% stack(3) is the main FieldTrip function that we are interested in
stack = stack(3);

% add information about the FieldTrip and MATLAB version used to the configuration
try
  cfg.callinfo.fieldtrip = ft_version();
catch
  cfg.callinfo.fieldtrip = 'unknown';
end
cfg.callinfo.matlab    = version();
cfg.callinfo.computer  = lower(computer); % for example maci64, glnx86, ...

% add information about the execution environment to the configuration
cfg.callinfo.hostname = gethostname();
cfg.callinfo.user     = getusername();
cfg.callinfo.pwd      = pwd;
cfg.callinfo.calltime = clock();

% add information about the function filename and revision to the configuration
cfg.version.name = stack.file;
clear stack

% the revision number is maintained by SVN in the revision variable in the calling function
if ~exist('revision', 'var')
  cfg.version.id   = 'unknown';
else
  cfg.version.id   = revision;
end

ftohDiW7th_FuncTimer = tic();
ftohDiW7th_FuncMem   = memtic();

