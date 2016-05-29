function spisop_hypvals(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, coreParameterFileName, parameterFileName)
% determine sleep scoring values for sleep table, assumes there is at least
% one valid epoch of S2 or S3 or S4 or REM
% Copyright Frederik D. Weber

spisop_hypvals_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, readParameterFileAsList(pathInputFolder,coreParameterFileName), readParameterFileAsList(pathInputFolder,parameterFileName));

end