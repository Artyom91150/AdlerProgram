function varargout = MAIN(varargin)
% MAIN MATLAB code for MAIN.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MAIN_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MAIN_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MAIN

% Last Modified by GUIDE v2.5 01-Nov-2020 10:35:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MAIN_OpeningFcn, ...
                   'gui_OutputFcn',  @MAIN_OutputFcn, ...
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


% --- Executes just before MAIN is made visible.
function MAIN_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MAIN (see VARARGIN)

% Choose default command line output for MAIN
handles.output = hObject;

handles.PlotParameters = InitializePlotParameters(handles);

handles.ArrayName = 'IC_Param_Array_0';
load('Map1');
handles.Map = Map;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MAIN wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MAIN_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

InitializeData(handles);
DrawMap(handles);

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Устанавливаются начальные значения элементов
function InitializeData(handles)
    
handles.slider1.Value = handles.Map.MapParam(3);
handles.slider1.Min = handles.Map.MapParam(3);
handles.slider1.Max = handles.Map.MapParam(4);

handles.slider2.Value = handles.Map.MapParam(5);
handles.slider2.Min = handles.Map.MapParam(5);
handles.slider2.Max = handles.Map.MapParam(6);

handles.edit2.String = handles.Map.MapParam(5);
handles.edit1.String = handles.Map.MapParam(3);

function PlotParameters = InitializePlotParameters(handles)
SystemParameters = struct('Gamma', [1.01, 1.01], 'Sigma', 0, 'd', 0);
FigureParameters = struct();
opts = odeset('RelTol',1e-9,'AbsTol',1e-12, 'MaxStep', 1e-1);
NumericalParameters = struct('InitialCondition', [0, 0], 'Tspan', str2double(handles.edit6.String), 'Options', opts);
ListOfComponents = struct('AreaOfD', 1, ...
                          'Isocline', handles.checkbox3.Value, ...
                          'Separatrix', handles.checkbox4.Value, ...
                          'MainTrajectories', handles.checkbox1.Value, ...
                          'SupportTrajectory', handles.checkbox2.Value, ...
                          'Equilibria', handles.checkbox5.Value, ...
                          'VectorField', handles.checkbox6.Value);

PlotParameters = struct('SystemParameters', SystemParameters, 'FigureParameters', FigureParameters, 'NumericalParameters', NumericalParameters, 'ListOfComponents', ListOfComponents);

function DrawMap(handles)

Visualization(handles.Map, handles.axes1);
%Sigma_d_Equations(handles.PlotParameters.SystemParameters.Gamma, handles.axes1);
plot(handles.axes1, handles.PlotParameters.SystemParameters.Sigma, handles.PlotParameters.SystemParameters.d, 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 0 0], 'LineWidth', 2);


function DrawPhaseSpace(ax, PlotParameters)
Phase_Space_Main(ax, PlotParameters);

function DrawTimeSeries(ax, PlotParameters)
Time_Series_Main(ax, PlotParameters);

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Value = get(hObject,'Value');

switch(Value)
   
    case 1
       load('Map1'); 
       delta = 0;
       %Массив для чтения и записи параметров и НУ
       %handles.ArrayName = 'IC_Param_Array_0_05';
    
    case 2
       load('Map_0_05_6');
       delta = 0.05;
       %Массив для чтения и записи параметров и НУ
       handles.ArrayName = 'IC_Param_Array_0_05_2';
        
    case 3
       load('Map_0_5_1'); 
       delta = 0.5;
       %Массив для чтения и записи параметров и НУ
       handles.ArrayName = 'IC_Param_Array_0_5_2';
    
end

handles.Map = Map;
handles.IC = [0, 0];
handles.PlotParameters.SystemParameters.Gamma = [1.01, 1.01 + delta];
guidata(hObject, handles);

InitializeData(handles);
DrawMap(handles);


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.SystemParameters.Sigma = str2double(get(hObject,'String'));
handles.slider1.Value = handles.PlotParameters.SystemParameters.Sigma;
DrawMap(handles);
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PlotParameters.SystemParameters.d = str2double(get(hObject,'String'));
handles.slider2.Value = handles.PlotParameters.SystemParameters.d;
DrawMap(handles);
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.SystemParameters.Sigma = get(hObject,'Value');
handles.edit1.String = num2str(handles.PlotParameters.SystemParameters.Sigma);
DrawMap(handles);
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.SystemParameters.d = get(hObject,'Value');
handles.edit2.String = num2str(handles.PlotParameters.SystemParameters.d);
DrawMap(handles);
guidata(hObject, handles);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
    ax = gca;
    fig = gcf;

Mode = [handles.checkbox1.Value, handles.checkbox2.Value];

if Mode(2) == 1
    %Параметры фигуры
    fig.Position = [592 482 1214 596];

    
    ax1 = subplot(1, 2, 1);
    axis(ax1, 'square');
    hold(ax1, 'on');
    DrawPhaseSpace(ax1, handles.PlotParameters);
    
    ax2 = subplot(1, 2, 2);
    axis(ax2, 'square');
    hold(ax2, 'on');
    DrawTimeSeries(ax2, handles.PlotParameters);
else
    fig.Position = [234 27 1026 963];
    axis(ax, 'square');
    hold(ax, 'on');
    DrawPhaseSpace(ax, handles.PlotParameters);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load(handles.ArrayName);

SpecialPlotParam = InitializePlotParameters(handles);

FontSize = 34;

%for i = 2:47
k = 1;
smth = ['(а)'; '(б)'; '(в)'; '(г)'];
for i = [3, 4, 7, 8]
    
    %handles.text9.String = [num2str(i-2), '/15'];

    if isempty(IC_Param_Array{1, i}) == 0
    
    
        Name = IC_Param_Array{1, i};
        SpecialPlotParam.SystemParameters.Sigma = IC_Param_Array{2, i};
        SpecialPlotParam.SystemParameters.d = IC_Param_Array{3, i};
        SpecialPlotParam.SystemParameters.Gamma = handles.PlotParameters.SystemParameters.Gamma;
        SpecialPlotParam.NumericalParameters.InitialCondition = [IC_Param_Array{4, i}, IC_Param_Array{5, i}];
        SpecialPlotParam.NumericalParameters.Tspan = IC_Param_Array{6, i};

        figure;
        ax = gca;
        fig = gcf;

        %if SpecialPlotParam.NumericalParameters.InitialCondition(1) ~= 0
            %Параметры фигуры
            %fig.Position = [1 41 1920 963];
            
            %fig.Position = [906 610 1214 612];
            fig.Position = [802 454 820 705];
            
            fig.FileName = Name;
            %fig.PaperPosition = [-5.3353 6.9654 31.6706 15.7692];
            %fig.PaperSize = [21.0000 29.7000];

%             ax1 = subplot(1, 2, 1);
%             axis(ax1, 'square');
%             hold(ax1, 'on');
%             DrawPhaseSpace(ax1, SpecialPlotParam);

            %ax2 = subplot(1, 2, 2);
            ax2 = gca;
            axis(ax2, 'square');
            hold(ax2, 'on');
            DrawTimeSeries(ax2, SpecialPlotParam)
            
            %XText1 = -ax1.XLim(2) / 5;
            XText2 = -ax2.XLim(2) / 5;
            YText = -1.5;
            
            text(ax2, XText2, -0.8, smth(k, :), 'FontSize', 40, 'interpreter', 'latex');
            k = k + 1;

%             if str2double(Name(length(Name))) == 1
%                 text(ax1, XText1, YText,'(а)', 'FontSize',FontSize, 'interpreter', 'latex');
%                 text(ax2, XText2, YText,'(б)', 'FontSize',FontSize, 'interpreter', 'latex');
%             else if str2double(Name(length(Name))) == 2
%                     text(ax1, XText1, YText,'(в)', 'FontSize',FontSize, 'interpreter', 'latex');
%                     text(ax2, XText2, YText,'(г)', 'FontSize',FontSize, 'interpreter', 'latex');
%                 else if str2double(Name(length(Name))) == 3
%                         text(ax1, XText1, YText,'(д)', 'FontSize',FontSize, 'interpreter', 'latex');
%                         text(ax2, XText2, YText,'(е)', 'FontSize',FontSize, 'interpreter', 'latex');
%                     else if str2double(Name(length(Name))) == 4
%                             text(ax1, XText1, YText,'(ж)', 'FontSize',FontSize, 'interpreter', 'latex');
%                             text(ax2, XText2, YText,'(з)', 'FontSize',FontSize, 'interpreter', 'latex');
%                         end
%                     end
%                 end
%             end
    
%         else
%         fig.Position = [234 27 1026 963];
%         axis(ax, 'square');
%         hold(ax, 'on');
%         DrawPhaseSpace(ax, SpecialPlotParam)
        %end
        
        %saveas(fig, Name, 'epsc2');
        %saveas(fig, Name, 'epsc');
        %saveas(fig, Name, 'svg');
        pause(0.1);
        %close(fig);

    end
    
    handles.text9.String = 'Готово';

end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load(handles.ArrayName);
i = 2;

while isempty(IC_Param_Array{1, i}) == 0
   i = i + 1; 
end

IC_Param_Array{1, i} = handles.edit3.String;
IC_Param_Array{2, i} = str2double(handles.edit1.String);
IC_Param_Array{3, i} = str2double(handles.edit2.String);
IC_Param_Array{6, i} = handles.PlotParameters.NumericalParameters.Tspan;

if handles.checkbox2.Value == 1
    IC_Param_Array{4, i} = str2double(num2str(handles.PlotParameters.NumericalParameters.InitialCondition(1)));
    IC_Param_Array{5, i} = str2double(num2str(handles.PlotParameters.NumericalParameters.InitialCondition(2)));
    
else
    IC_Param_Array{4, i} = 0;
    IC_Param_Array{5, i} = 0;
end


save(handles.ArrayName, 'IC_Param_Array');


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.ListOfComponents.MainTrajectories = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.ListOfComponents.SupportTrajectory = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.ListOfComponents.Isocline = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.ListOfComponents.Separatrix = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.ListOfComponents.Equilibria = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    
handles.PlotParameters.NumericalParameters.InitialCondition(1) = str2double(get(hObject,'String'));
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PlotParameters.NumericalParameters.InitialCondition(2) = str2double(get(hObject,'String'));
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Mode = [handles.checkbox1.Value, handles.checkbox2.Value];

if Mode(2) == 1
    
    if handles.radiobutton1.Value == 1
        handles.PlotParameters.NumericalParameters.InitialCondition = ginput;
        handles.edit4.String = num2str(handles.PlotParameters.NumericalParameters.InitialCondition(1));
        handles.edit5.String = num2str(handles.PlotParameters.NumericalParameters.InitialCondition(2));
    else
        handles.PlotParameters.NumericalParameters.InitialCondition = [str2double(handles.edit4.String), str2double(handles.edit5.String)];
    end
    guidata(hObject, handles);
    
    DrawPhaseSpace(handles.axes2, handles.PlotParameters);
    DrawTimeSeries(handles.axes3, handles.PlotParameters);
else
    DrawPhaseSpace(handles.axes2, handles.PlotParameters);
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PlotParameters.ListOfComponents.VectorField = get(hObject,'Value');
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox6



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.PlotParameters.NumericalParameters.Tspan = str2double(get(hObject,'String'));
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
