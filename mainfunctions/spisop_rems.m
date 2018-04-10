function spisop_rems(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, coreParameterFileName, parameterFileName)
% rems and related phenomena detection
% Copyright Frederik D. Weber

spisop_rems_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, readParameterFileAsList(pathInputFolder,coreParameterFileName), readParameterFileAsList(pathInputFolder,parameterFileName));

end