function listOfParameters = readParameterFileAsList(pathInputFolder,parameterFileName)

if exist([pathInputFolder filesep parameterFileName],'file') ~= 2
    error(['Function Parameter file ' [pathInputFolder filesep parameterFileName] ' does not exist. Check for correct path and if file exists in it.'])
end

listOfParameters = read_mixed_csv([pathInputFolder filesep parameterFileName],',');

end