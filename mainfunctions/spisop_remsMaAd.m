function spisop_remsMaAd(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, coreParameterFileName, parameterFileName)
% rems and related phenomena detection
% Copyright Frederik D. Weber

spisop_remsMaAd_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, readParameterFileAsList(pathInputFolder,coreParameterFileName), readParameterFileAsList(pathInputFolder,parameterFileName));

end