function spisop_eventCooccurrence(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, parameterFileName)
% discover co-occurrence of test and target events, i.e. if test events
% fall within a defined timewindow arround target events.
% Copyright Frederik D. Weber

spisop_eventCooccurrence_l1(pathInputFolder, pathOutputFolder, ouputFilesPrefixString, readParameterFileAsList(pathInputFolder,parameterFileName));

end