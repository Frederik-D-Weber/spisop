%  The function should take the cfg-structure as input and should give a
%  NxM matrix with M equal to or larger than 3) in the same format as
%  "trl" as the output. You can add extra custom fields to the
%  configuration structure to pass as arguments to your own trialfun.
%  Furthermore, inside the trialfun you can use the FT_READ_EVENT
%  function to get the event information from your data file.
%
%   The trial definition "trl" is an Nx3 matrix, N is the number of trials.
%   The 
% first column contains the sample-indices of the begin of each trial relative to the begin of the raw data, 
% the second column contains the sample-indices of the end of each trial,
% the third column contains the offset of the trigger with respect to the trial. 
%     An offset of 0  means that the first sample of the trial corresponds to the trigger.
%     A positive offset indicates that the first sample is later than the trigger, 
%     a negative offset indicates that the trial begins before the trigger.
%  
%   The trial definition "trl" can contain additional columns besides the
%   required three that represend begin, end and offset. These additional
%   columns can be used by a custom trialfun to provide numeric information
%   about each trial such as trigger codes, response latencies, trial
%   type and response correctness. The additional columns of the "trl"
%   matrix will be represented in data.trialinfo after FT_PREPROCESSING.

function trl = trialfun_spd_ROIs(cfg)
  trl = [cfg.roiBegins cfg.roiEnds zeros(length(cfg.roiEnds),1)];
end