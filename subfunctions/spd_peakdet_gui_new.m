function varargout = spd_peakdet_gui_new(varargin)
% SPD_PEAKDET_GUI_NEW MATLAB code for spd_peakdet_gui_new.fig
%      SPD_PEAKDET_GUI_NEW, by itself, creates a new SPD_PEAKDET_GUI_NEW or raises the existing
%      singleton*.
%
%      H = SPD_PEAKDET_GUI_NEW returns the handle to a new SPD_PEAKDET_GUI_NEW or the handle to
%      the existing singleton*.
%
%      SPD_PEAKDET_GUI_NEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPD_PEAKDET_GUI_NEW.M with the given input arguments.
%
%      SPD_PEAKDET_GUI_NEW('Property','Value',...) creates a new SPD_PEAKDET_GUI_NEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spd_peakdet_gui_new_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spd_peakdet_gui_new_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spd_peakdet_gui_new

% Last Modified by GUIDE v2.5 18-Jun-2015 12:46:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spd_peakdet_gui_new_OpeningFcn, ...
                   'gui_OutputFcn',  @spd_peakdet_gui_new_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before spd_peakdet_gui_new is made visible.
function spd_peakdet_gui_new_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spd_peakdet_gui_new (see VARARGIN)

% Choose default command line output for spd_peakdet_gui_new
handles.output = varargin{3};

% UIWAIT makes spd_peakdet_gui_new wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.pFreq = varargin{1};
handles.pFreqDet = varargin{1};
handles.pPower = varargin{2};
handles.pPowerDet = varargin{2};
handles.pPowerChan = varargin{3};
handles.pPowerChanDet = varargin{3};
handles.pPowerChanLabels = varargin{4};
%handles.pPower = log10(varargin{2});
handles.freqPeaks = varargin{5};
handles.freqPeaksDet = varargin{5};
handles.powerPeaks = varargin{6};
handles.powerPeaksDet = varargin{6};
%handles.powerPeaks = log10(varargin{4});
%handles.powerPeaksDet = log10(varargin{4});
handles.tempiDataString = varargin{7};
handles.nPeaks = length(handles.freqPeaks);
handles.numberOfPossiblePeaks = 2;
handles.ouputFilesPrefixString = varargin{8};
handles.datasetnum = varargin{9};

for k = 1:handles.nPeaks
set(handles.(['edit_peak_' num2str(k)]),'String',num2str(handles.freqPeaks(k)));
end
if handles.numberOfPossiblePeaks > handles.nPeaks
    for k = (handles.nPeaks+1):handles.numberOfPossiblePeaks
    set(handles.(['edit_peak_' num2str(k)]),'String','UNDEF');
    end
end

handles.tempColors = hsv(length(handles.freqPeaks));

handles.axes_freq_power = plot(handles.pFreq,handles.pPower,'Color',[0 0 0],'LineWidth',1.5,'DisplayName','mean(Chan)');
    hold on
    temp_colors = hsv(length(handles.pPowerChanLabels));
    for iChan = 1:length(handles.pPowerChanLabels) 
        temp_plot = plot(handles.pFreq,handles.pPowerChan(iChan,:),'Color',temp_colors(iChan,:),'LineWidth',0.5,'DisplayName',char(strrep(handles.pPowerChanLabels(iChan),'_','\_')));
        
    end
    legend(gca,'show');
    legend('boxoff');
    for iENP = 1:length(handles.freqPeaks)
        plot(handles.freqPeaks(iENP),handles.powerPeaks(iENP),'o','MarkerFaceColor',handles.tempColors(iENP,:))
        %line([handles.freqPeaks(iENP) handles.freqPeaks(iENP)],[0 handles.powerPeaks(iENP)],'LineStyle','--','Color',[0 0 0])
        vline(handles.freqPeaks(iENP),'LineStyle','--','Color',[0 0 0])
    end
    title([' power spectrum' ' dataset ' handles.tempiDataString]);
    xlabel('Frequency [Hz]');
    ylabel('Power [P^2, e.g. µV^2]');
    hold off
    
% Update handles structure
guidata(hObject, handles);

uiwait(handles.figure1); 

%set(handles.axes_freq_power,'XDataSource',pFreq,'YDataSource',pPower);


% --- Outputs from this function are returned to the command line.
function varargout = spd_peakdet_gui_new_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.freqPeaks;
delete(handles.figure1);

% --- Executes on button press in pushbutton_accept.
function pushbutton_accept_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(gcf, 'PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
print([handles.ouputFilesPrefixString 'power_spectrum_dataset_' num2str(handles.datasetnum)],'-dpng','-r300');
uiresume(handles.figure1);




function edit_peak_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_peak_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_peak_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_peak_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_peak_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_peak_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_peak_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_peak_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_peak_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_peak_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_peak_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_peak_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_refresh.
function pushbutton_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for k = 1:handles.nPeaks
    handles.freqPeaks(k) = str2double(get(handles.(['edit_peak_' num2str(k)]),'String'));
end
handles.freqPeaks = sort(handles.freqPeaks,'ascend');
for k = 1:handles.nPeaks
    temp = handles.pPower(find(handles.pFreq >= handles.freqPeaks(k),1,'first'));
    if length(temp) < 1
         handles.powerPeaks(k) = max(handles.pFreq);
    else
        handles.powerPeaks(k) = temp;
    end
end


handles.axes_freq_power = plot(handles.pFreq,handles.pPower,'Color',[0 0 0],'LineWidth',1.5,'DisplayName','mean(Chan)');
    hold on
    temp_colors = hsv(length(handles.pPowerChanLabels));
    for iChan = 1:length(handles.pPowerChanLabels) 
        temp_plot = plot(handles.pFreq,handles.pPowerChan(iChan,:),'Color',temp_colors(iChan,:),'LineWidth',0.5,'DisplayName',char(strrep(handles.pPowerChanLabels(iChan),'_','\_')));  
    end
    legend(gca,'show');
    legend('boxoff');
    for iENP = 1:length(handles.freqPeaks)
        if (handles.freqPeaksDet(iENP)>=min(handles.pFreq)) && (handles.freqPeaksDet(iENP)<=max(handles.pFreq))
            plot(handles.freqPeaksDet(iENP),handles.powerPeaksDet(iENP),'r*')
        end
        if (handles.freqPeaks(iENP)>=min(handles.pFreq)) && (handles.freqPeaks(iENP)<=max(handles.pFreq))
            plot(handles.freqPeaks(iENP),handles.powerPeaks(iENP),'o','MarkerFaceColor',handles.tempColors(iENP,:))
            %line([handles.freqPeaks(iENP) handles.freqPeaks(iENP)],[0 handles.powerPeaks(iENP)],'LineStyle','--','Color',[0 0 0])
            vline(handles.freqPeaks(iENP),'LineStyle','--','Color',[0 0 0])
        end
    end
    title([' power spectrum' ' dataset ' handles.tempiDataString]);
    xlabel('Frequency (Hz)');
    ylabel('Power [P^2, e.g. µV^2]');
    hold off
    
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.pFreq = handles.pFreqDet;
handles.pPower = handles.pPowerDet;
handles.pPowerChan = handles.pPowerChanDet;

for k = 1:handles.nPeaks
    handles.freqPeaks(k) = handles.freqPeaksDet(k);
    handles.powerPeaks(k) = handles.pPower(find(handles.pFreq >= handles.freqPeaks(k),1,'first'));
end

handles.axes_freq_power = plot(handles.pFreq,handles.pPower,'Color',[0 0 0],'LineWidth',1.5,'DisplayName','mean(Chan)');
    hold on
    temp_colors = hsv(length(handles.pPowerChanLabels));
    for iChan = 1:length(handles.pPowerChanLabels) 
        temp_plot = plot(handles.pFreq,handles.pPowerChan(iChan,:),'Color',temp_colors(iChan,:),'LineWidth',0.5,'DisplayName',char(strrep(handles.pPowerChanLabels(iChan),'_','\_')));
    end
    legend(gca,'show');
    legend('boxoff');
    for iENP = 1:length(handles.freqPeaks)
        plot(handles.freqPeaksDet(iENP),handles.powerPeaksDet(iENP),'r*')
        plot(handles.freqPeaks(iENP),handles.powerPeaks(iENP),'o','MarkerFaceColor',handles.tempColors(iENP,:))
        %line([handles.freqPeaks(iENP) handles.freqPeaks(iENP)],[0 handles.powerPeaks(iENP)],'LineStyle','--','Color',[0 0 0])
        vline(handles.freqPeaks(iENP),'LineStyle','--','Color',[0 0 0])
    end
    title([' power spectrum' ' dataset ' handles.tempiDataString]);
    xlabel('Frequency (Hz)');
    ylabel('Power [P^2, e.g. µV^2]');
    hold off
    
for k = 1:handles.nPeaks
    set(handles.(['edit_peak_' num2str(k)]),'String',num2str(handles.freqPeaks(k)));
end



% Update handles structure
guidata(hObject, handles);

% % --- Executes when user attempts to close figure1.
% function figure1_CloseRequestFcn (hObject, eventdata, handles)
% % hObject handle to figure1 (see GCBO)
% % eventdata reserved - to be defined in a future version of MATLAB
% % handles structure with handles and user data (see GUIDATA)
% % Hint: delete(hObject) closes the figure
% uiresume(handles.figure1);
% % delete(hObject); wurde hier entfernt!

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
% delete(hObject); wurde hier entfernt!


% --- Executes on button press in pushbutton_skip.
function pushbutton_skip_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_skip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for k = 1:handles.nPeaks
    handles.freqPeaks(k) = -1;
end
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_zoomFromLeft.
function pushbutton_zoomFromLeft_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_zoomFromLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.pFreq = handles.pFreq(round(length(handles.pFreq)*0.05+1):end);
handles.pPower = handles.pPower(round(length(handles.pPower)*0.05+1):end);
handles.pPowerChan = handles.pPowerChan(:,round(size(handles.pPowerChan,2)*0.05+1):end);

pushbutton_refresh_Callback(hObject, eventdata, handles)




% --- Executes on button press in pushbutton_zoomFromRight.
function pushbutton_zoomFromRight_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_zoomFromRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.pFreq = handles.pFreq(1:end-round(length(handles.pFreq)*0.05+1));
handles.pPower = handles.pPower(1:end-round(length(handles.pPower)*0.05+1));
handles.pPowerChan = handles.pPowerChan(:,1:end-round(size(handles.pPowerChan,2)*0.05+1));

pushbutton_refresh_Callback(hObject, eventdata, handles)
