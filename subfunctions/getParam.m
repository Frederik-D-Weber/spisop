function str = getParam(parameterstring,listOfParameters)
    str = listOfParameters{find(strcmp(listOfParameters(:,1),parameterstring),1,'first'),2};
end