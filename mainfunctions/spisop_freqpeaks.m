function spisop_freqpeaks(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, coreParameterFileName, parameterFileName)
% determine peaks in frequency band of e.g. spindles or SO etc. (depends on parameters)
% Copyright Frederik D. Weber

spisop_freqpeaks_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, readParameterFileAsList(pathInputFolder,coreParameterFileName), readParameterFileAsList(pathInputFolder,parameterFileName));

end
