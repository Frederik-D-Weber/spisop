function [ iDatas, doAsInParameterFile, useDummyDataset] = dataSetInputDialog(listOfDatasetsPaths,iDatas,selectionType,initDataIndex)

if initDataIndex > numel(listOfDatasetsPaths)
    initDataIndex = numel(listOfDatasetsPaths);
end

initDataIndex = initDataIndex + 1;

doAsInParameterFile = false;
temp_has_decided = false;
temp_selection_items = listOfDatasetsPaths;
useDummyDataset = false;
if strcmp(selectionType,'single')
    temp_selection_items = ['dummy' ; temp_selection_items ];
end
set(0,'units','pixels')
%Obtains this pixel information
temp_screensize = get(0,'screensize');
while ~temp_has_decided
    temp_listSize = floor(temp_screensize(3:4)/2);
    promt = 'Select muliple datasets';
    if strcmp(selectionType,'single')
        promt = 'Select a dataset';    end
    [dataset_selection,dataset_selection_ok] = listdlg('PromptString',promt,...
        'SelectionMode',selectionType,...
        'ListString',temp_selection_items,...
        'InitialValue',initDataIndex,...
        'ListSize',temp_listSize,...
        'CancelString','Open all as stated in Parameterfile');
    if dataset_selection_ok
        if strcmp(selectionType,'single')
            iDatas = dataset_selection(1)-1;
        else
            iDatas = dataset_selection;
        end
        temp_has_decided = true;
    else
        % Construct a questdlg with three options
        choice = questdlg('Really do NOT WANT to select one dataset and ignore parameterfile definitions?', ...
            'How to proceed?', ...
            'Like defined in Parameter file','I changed my mind, go back','I changed my mind, go back');
        % Handle response
        switch choice
            case 'I changed my mind, go back'
                %temp_has_decided = false;
            otherwise
                temp_has_decided = true;
                doAsInParameterFile = true;
        end
    end
end
if all(iDatas == 0)
    useDummyDataset = true;
end
end

