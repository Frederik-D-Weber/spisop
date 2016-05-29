function spisop_nonEvents(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, coreParameterFileName, parameterFileName)
% discover correspondent non-events according to new events and hypnogram
% Copyright Frederik D. Weber

spisop_nonEvents_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, readParameterFileAsList(pathInputFolder,coreParameterFileName), readParameterFileAsList(pathInputFolder,parameterFileName));

end


