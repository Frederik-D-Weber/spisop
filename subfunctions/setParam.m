function lop = setParam(parameterstring,listOfParameters,strParameter)
    lop = listOfParameters;
    lop{find(strcmp(listOfParameters(:,1),parameterstring),1,'first'),2} = strParameter;
end