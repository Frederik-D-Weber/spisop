function str = getParam(parameterstring,listOfParameters)
ind = find(strcmp(listOfParameters(:,1),parameterstring),1,'first');
if isempty(ind)
        error(['The parameter string ' parameterstring ' is was not found in the parameter file, please check\n Parameter file list dimension is ' num2str(size([listOfParameters]))])
end
    str = listOfParameters{ind,2};
end