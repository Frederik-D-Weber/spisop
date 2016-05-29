function res = pretestHeaderForPersistentSampleFrequencies(IgnoreDataSetHeader,iDatas,listOfDatasetHeaderPaths,listOfChannelsOfInterest,FrqOfSmplWished)
if strcmp(IgnoreDataSetHeader,'no')
    fprintf('Pre-test if header sample frequencies in datasets are lower than to wished frequency\n');
    hasNotDataRequiredFrqOfSmplWished = logical(zeros(1,length(iDatas)));
    headerFreqs =  zeros(1,length(iDatas));
    hasNoChannels = logical(zeros(1,length(iDatas)));
    tempConseciDatas = 1:length(iDatas);
    for tempConseciData = tempConseciDatas
    %parfor tempConseciData = tempConseciDatas

        iData = iDatas(tempConseciData);
        headerPath = listOfDatasetHeaderPaths{iData};
        hdr = [];
        hdr = ft_read_header(headerPath);
        if (FrqOfSmplWished > hdr.Fs)
            hasNotDataRequiredFrqOfSmplWished(iData) = true;
        end
        headerFreqs(iData) = hdr.Fs;
        
        tempChannelsOfInterest = listOfChannelsOfInterest(iData,:);
        tempChannelsOfInterest = tempChannelsOfInterest(~(cellfun(@isempty,tempChannelsOfInterest)));
        tempChannels = ft_channelselection(tempChannelsOfInterest, hdr.label);
        
        if (isempty(tempChannels))
            hasNoChannels(iData) = true;
        end

    end
    if any(hasNotDataRequiredFrqOfSmplWished)
        tempdatas = 1:length(iDatas);
        error(['datasets ' num2str(tempdatas(hasNotDataRequiredFrqOfSmplWished)) ': designated wished frequency of ' num2str(FrqOfSmplWished) ' not supported by respective data frequencies of ' num2str(headerFreqs(hasNotDataRequiredFrqOfSmplWished)) ' Hz, choose lower FrqOfSmplWished that is supported by all those datasets !']);
    end
    if any(hasNoChannels)
        tempdatas = 1:length(iDatas);
        error(['datasets ' num2str(tempdatas(hasNoChannels)) ': has not even one channel selected, choose channels of interest new, or exclude these datasets from the analysis!']);
    end
end
res = true;
end