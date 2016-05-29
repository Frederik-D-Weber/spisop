function spisop_manipulateParam(pathInputFolder, inParamFilename, outParamFilename, param, value_str)



if exist([pathInputFolder filesep inParamFilename],'file') ~= 2
    error(['input Parameter File ' [pathInputFolder filesep inParamFilename] ' does not exist. Check for correct path and if file exists in it.'])
end

temp_listOfinParameters = read_mixed_csv([pathInputFolder filesep inParamFilename],',','includecomments');

try
    getParam(param,temp_listOfinParameters);
catch err
    error(['the Parameter ' param ' was not found in the input Parameter File ' [pathInputFolder filesep inParamFilename]]);
end
temp_listOfOutParameters = setParam(param,temp_listOfinParameters,value_str);

writetable(cell2table(temp_listOfOutParameters),[pathInputFolder filesep outParamFilename],'Delimiter',',','FileType','text','WriteVariableNames',0,'WriteRowNames',0);

end
