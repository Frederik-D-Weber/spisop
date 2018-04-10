function spisop_browser(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, coreParameterFileName, parameterFileName)
% a browser for data an events
% Copyright Frederik D. Weber

pathInputFolder_coreparam = pathInputFolder;
pathInputFolder_param = pathInputFolder;

[temppath_coreparam, tempname_coreparam, tempext_coreparam] = fileparts(coreParameterFileName);
if ~strcmp(temppath_coreparam,'') && strcmp(pathInputFolder,'')
    pathInputFolder_coreparam = temppath_coreparam;
    coreParameterFileName = [tempname_coreparam, tempext_coreparam];
end
[temppath_param, tempname_param, tempext_param] = fileparts(parameterFileName);

if ~strcmp(temppath_param,'') && strcmp(pathInputFolder,'')
    pathInputFolder_param = temppath_param;
    parameterFileName = [tempname_param, tempext_param];
end

listOfCoreParameters = readParameterFileAsList(pathInputFolder_coreparam,coreParameterFileName);
listOfParameters = readParameterFileAsList(pathInputFolder_param,parameterFileName);


spisop_browser_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, listOfCoreParameters, listOfParameters);

end