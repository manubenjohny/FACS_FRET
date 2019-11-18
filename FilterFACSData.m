function FilterFACSData(varargin)
% [OutData ROIStruct] = FilterFACSData(0, FileName) or simpley FilterFACSData for mainFxn
% [OutData ROIStruct] = FilterFACSData(1, 'CallBackFxn',variables) for callback
if isempty(varargin)
    varargin{1} = 0;
end
if varargin{1}==0
    InitializeGUI(varargin{2:end});
else
    if ischar(varargin{2}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
        try
            if (nargout)
                [varargout{1:nargout}] = feval(varargin{2}); % FEVAL switchyard
            else
                feval(varargin{2:end}); % FEVAL switchyard
            end
        catch EM
            rethrow(EM);
        end
    end
end

function InitializeGUI(varargin)

if nargin == 0
    [fname pathname] = uigetfile('*.fcs', 'Open Raw Data');  
    try
        FACSData = ReadFACSFile([pathname fname]);
    catch EM
        rethrow(EM);
        return;
    end
    GUIDAT.WorkingPath = pathname;
    GUIDAT.Filename = fname;
    RawData = FACSData.Data;
    Header = FACSData.Header;
    bFN = 1;
elseif nargin == 1    
    % initialize variables etc;
    FACSData = varargin{1};
    RawData = FACSData.Data;
    Header = FACSData.Header;
else 
    OutData = [];

    return;   
end

% Setup Buttons Etc
if bFN>0
    MyHandles.FG = figure('name',['Filter Raw FCS Data | ' pathname fname], 'NumberTitle','off','units','normalized', 'position', [0.1 0.1 0.8 0.4],'ToolBar','none','MenuBar','none','color', [1 1 1],'KeyPressFcn',@keypresser);
else
    MyHandles.FG = figure('name','Filter Raw FCS Data', 'NumberTitle','off','units','normalized', 'position', [0.1 0.1 0.8 0.4],'ToolBar','none','MenuBar','none','color', [1 1 1],'KeyPressFcn',@keypresser);
end
MyHandles.Ax(1) = axes('units','normalized','position',[0.05 0.12 0.25 0.7]);
MyHandles.Ax(2) = axes('units','normalized','position',[0.35 0.12 0.25 0.7]);
MyHandles.Ax(3) = axes('units','normalized','position',[0.65 0.12 0.25 0.7]);

for j = 1:3
    MyHandles.PBR(j) = uicontrol('style','pushbutton','String', 'Pick ROI','units','normalized','position',[0.05+0.3*(j-1) 0.9 0.25 0.05],'Callback',['FilterFACSData(1,''PickROI'',' num2str(j) ')']);
end

MyHandles.ChBR1 = uicontrol('style','checkbox','String','In R1','units','normalized', 'position', [0.55 0.85 0.05 0.05],'Callback',['FilterFACSData(1,''RePlot'')']);
MyHandles.ChBR2 = uicontrol('style','checkbox','String','In R1','units','normalized', 'position', [0.8 0.85 0.05 0.05],'Callback',['FilterFACSData(1,''RePlot'')']);
MyHandles.ChBR3 = uicontrol('style','checkbox','String','In R2','units','normalized', 'position', [0.85 0.85 0.05 0.05],'Callback',['FilterFACSData(1,''RePlot'')']);
MyHandles.PBPlotF = uicontrol('style','pushbutton','String', 'PlotF','units','normalized','position',[0.9 0.9 0.05 0.05],'Callback',['FilterFACSData(1,''PlotF'')']);
MyHandles.PBSave = uicontrol('style','pushbutton','String', 'Save','units','normalized','position',[0.9 0.85 0.05 0.05],'Callback',['FilterFACSData(1,''Save'')']);
MyHandles.PBSaveROI = uicontrol('style','pushbutton','String', 'Save ROI','units','normalized','position',[0.9 0.8 0.05 0.05],'Callback',['FilterFACSData(1,''SaveROI'')']);
MyHandles.PBLoadROI = uicontrol('style','pushbutton','String', 'Load ROI','units','normalized','position',[0.9 0.75 0.05 0.05],'Callback',['FilterFACSData(1,''LoadROI'')']);
MyHandles.PBAdvance = uicontrol('style','pushbutton','String', '>','units','normalized','position',[0.925 0.7 0.025 0.05],'Callback',['FilterFACSData(1,''AdvanceFile'')']);
MyHandles.PBRetreat = uicontrol('style','pushbutton','String', '<','units','normalized','position',[0.9 0.7 0.025 0.05],'Callback',['FilterFACSData(1,''RetreatFile'')']);
MyHandles.ProcessAll = uicontrol('style','pushbutton','String', 'ProcessAll','units','normalized','position',[0.9 0.65 0.05 0.05],'Callback',['FilterFACSData(1,''ProcessAllFiles'')']);


%Initialize Data Structure -- DATA has RawData and ROI
DATA.RawData = RawData; 
DATA.Header = Header;

DATA.ROI(1).Edge.X = [0];
DATA.ROI(1).Edge.Y = [0];
DATA.ROI(2).Edge.X = [0];
DATA.ROI(2).Edge.Y = [0];
DATA.ROI(3).Edge.X = [0];
DATA.ROI(3).Edge.Y = [0];
DATA.MASK = ones(length(RawData.FSCH),3);

FSCH_LIM = [0 min(trimmean(RawData.FSCH,99)+std(RawData.FSCH)*2.5,2.5E5)];
SSCH_LIM = [0 min(trimmean(RawData.SSCH,99)+std(RawData.SSCH)*2.5,2.5E5)];
FSCA_LIM = [0 min(trimmean(RawData.FSCA,99)+std(RawData.FSCA)*2.5,2.5E5)];
SSCA_LIM = [0 min(trimmean(RawData.SSCA,99)+std(RawData.SSCA)*2.5,2.5E5)];

axes(MyHandles.Ax(1));
MyHandles.DAT(1) = PlotDATA2D(RawData.FSCH,RawData.SSCH);
MyHandles.ROI(1) = line(0,0);
MyHandles.ROI(1) = PlotROI(DATA.ROI(1).Edge.X,DATA.ROI(1).Edge.Y,MyHandles.ROI(1));
set(MyHandles.Ax(1),'TickDir','out','XLim',FSCH_LIM,'YLim',SSCH_LIM)
xlabel('FSC-H');
ylabel('SSC-H');


axes(MyHandles.Ax(2));
MyHandles.DAT(2) = PlotDATA2D(RawData.FSCA,RawData.FSCH);
set(MyHandles.Ax(2),'TickDir','out','XLim',FSCA_LIM,'YLim',FSCH_LIM)
MyHandles.ROI(2) = line(0,0);
MyHandles.ROI(2) = PlotROI(DATA.ROI(2).Edge.X,DATA.ROI(2).Edge.Y,MyHandles.ROI(2));
xlabel('FSC-A');
ylabel('FSC-H');

axes(MyHandles.Ax(3));
MyHandles.DAT(3) = PlotDATA2D(RawData.SSCA,RawData.SSCH);
set(MyHandles.Ax(3),'TickDir','out','XLim',SSCA_LIM,'YLim',SSCH_LIM)
MyHandles.ROI(3) = line(0,0);
MyHandles.ROI(3) = PlotROI(DATA.ROI(3).Edge.X,DATA.ROI(3).Edge.Y,MyHandles.ROI(3));
xlabel('SSC-A');
ylabel('SSC-H');

GUIDAT.DATA = DATA;
GUIDAT.MyHandles = MyHandles;

set(MyHandles.FG ,'UserData',GUIDAT)



function PickROI(bVal)
GUIDAT = get(gcf,'UserData');
% draw ROI
axes(GUIDAT.MyHandles.Ax(bVal))
drawNewROI(bVal);
% set(gcf,'UserData',GUIDAT);
FilterData()
set(GUIDAT.MyHandles.ChBR1,'value',1);
RePlot()

function Handles = PlotDATA2D(XDATA,YDATA,Ax)
Handles.PTS = line('XData',XDATA,'YData',YDATA,'LineStyle','none','Marker','o','MarkerSize',0.5,'MarkerFaceColor',[0.4 0.4 0.4],'MarkerEdgeColor','none');
[n,c] = hist3([YDATA,XDATA], [100 100]*3);
hold on;
[~, Handles.CTR] = contour(c{2}, c{1}, n,20);
hold off;

function [MyL PolyEdge] = drawNewROI(bVal)
GUIDAT = get(gcf,'UserData');
MyL = GUIDAT.MyHandles.ROI(bVal);
j = 1;
bstop = 0;
% x = [0];
% y = [0];  
% MyL = line(x,y,'LineStyle','-','Color',[1 0 0],'LineWidth',3,'MarkerFaceColor',[1 0 0]); 
while ~bstop
    [xtemp,ytemp] = ginput(1);
    if ~isempty(xtemp)
        x(j) = xtemp; y(j) = ytemp;
    else
        bstop = 1;
    end
    set(MyL,'XData', x,'YData',y,'LineStyle','-','Color',[1 0 0],'Marker','o');
j = j+1;
end
x = [x, x(1)];
y = [y, y(1)];
set(MyL,'XData', x,'YData',y,'LineStyle','-','Color',[0 1 0],'Marker','o');
GUIDAT.DATA.ROI(bVal).Edge.X = x;
GUIDAT.DATA.ROI(bVal).Edge.Y = y;
set(gcf,'UserData',GUIDAT);


function handle = PlotROI(XData, YData,handle)
set(handle,'XData',XData,'YData', YData,'LineStyle','-','Color',[0 1 0],'LineWidth',3,'MarkerFaceColor',[0 1 0]); 


function FilterData()
GUIDAT = get(gcf,'UserData');
if length(GUIDAT.DATA.ROI(1).Edge.X)>2
    [IN ON] = inpolygon(GUIDAT.DATA.RawData.FSCH,GUIDAT.DATA.RawData.SSCH,GUIDAT.DATA.ROI(1).Edge.X,GUIDAT.DATA.ROI(1).Edge.Y);
    GUIDAT.DATA.MASK(:,1) = (IN | ON);
end
if length(GUIDAT.DATA.ROI(2).Edge.X)>2
    [IN ON] = inpolygon(GUIDAT.DATA.RawData.FSCA,GUIDAT.DATA.RawData.FSCH,GUIDAT.DATA.ROI(2).Edge.X,GUIDAT.DATA.ROI(2).Edge.Y);
    GUIDAT.DATA.MASK(:,2) = (IN | ON);
end
if length(GUIDAT.DATA.ROI(3).Edge.X)>2
    [IN ON] = inpolygon(GUIDAT.DATA.RawData.SSCA,GUIDAT.DATA.RawData.SSCH,GUIDAT.DATA.ROI(3).Edge.X,GUIDAT.DATA.ROI(3).Edge.Y);
    GUIDAT.DATA.MASK(:,3) = (IN | ON);
end
set(gcf,'UserData',GUIDAT);

function RePlot()
GUIDAT = get(gcf,'UserData');
MyHandles = GUIDAT.MyHandles;
DATA = GUIDAT.DATA;
if isfield(GUIDAT.DATA,'MASK')
    MASK = GUIDAT.DATA.MASK>0;

    axes(MyHandles.Ax(1));
    MyHandles.ROI(1) = PlotROI(DATA.ROI(1).Edge.X,DATA.ROI(1).Edge.Y,MyHandles.ROI(1));
    set(MyHandles.DAT(1).PTS,'XData',DATA.RawData.FSCH,'YData',DATA.RawData.SSCH);

    axes(MyHandles.Ax(2));
    if get(GUIDAT.MyHandles.ChBR1,'value')>0
       set(MyHandles.DAT(2).PTS,'XData',DATA.RawData.FSCA(MASK(:,1)),'YData',DATA.RawData.FSCH(MASK(:,1)));
       hold on;
       [n,c] = hist3([DATA.RawData.FSCH(MASK(:,1)),DATA.RawData.FSCA(MASK(:,1))], [50 50]);
       delete(MyHandles.DAT(2).CTR);
       [~, GUIDAT.MyHandles.DAT(2).CTR] = contour(c{2}, c{1}, n, 5);
       hold off;  
    else    
       set(MyHandles.DAT(2).PTS,'XData',DATA.RawData.FSCA,'YData',DATA.RawData.FSCH);
       hold on;
       [n,c] = hist3([DATA.RawData.FSCH,DATA.RawData.FSCA], [50 50]);
       delete(MyHandles.DAT(2).CTR);
       [~, GUIDAT.MyHandles.DAT(2).CTR] = contour(c{2}, c{1}, n, 5);
       hold off;  


    end
    MyHandles.ROI(2) = PlotROI(DATA.ROI(2).Edge.X,DATA.ROI(2).Edge.Y,MyHandles.ROI(2));


    axes(MyHandles.Ax(3));
    bM1 = get(GUIDAT.MyHandles.ChBR2,'value')>0;
    bM2 = get(GUIDAT.MyHandles.ChBR3,'value')>0;
    if (bM1 && bM2)
       TempMask = sum(MASK(:,1:2),2)>1; 
    elseif (bM1 && ~bM2)
       TempMask = MASK(:,1)>0; 
    elseif (~bM1 && bM2)
       TempMask = MASK(:,2)>0; 
    else   
       TempMask = ones(size(MASK(:,1)))>0; 
    end
    set(MyHandles.DAT(3).PTS,'XData',DATA.RawData.SSCA(TempMask),'YData',DATA.RawData.SSCH(TempMask));

       hold on;
       [n,c] = hist3([DATA.RawData.SSCH(TempMask),DATA.RawData.SSCA(TempMask)], [100 100]*3);
       delete(MyHandles.DAT(3).CTR);
       [~, GUIDAT.MyHandles.DAT(3).CTR] = contour(c{2}, c{1}, n, 5);
       hold off;  

    MyHandles.ROI(3) = PlotROI(DATA.ROI(3).Edge.X,DATA.ROI(3).Edge.Y,MyHandles.ROI(3));
    set(gcf,'UserData',GUIDAT)

else
FSCH_LIM = [0 min(trimmean(DATA.RawData.FSCH,99)+std(DATA.RawData.FSCH)*2,2.5E5)];
SSCH_LIM = [0 min(trimmean(DATA.RawData.SSCH,99)+std(DATA.RawData.SSCH)*2,2.5E5)];
FSCA_LIM = [0 min(trimmean(DATA.RawData.FSCA,99)+std(DATA.RawData.FSCA)*2,2.5E5)];
SSCA_LIM = [0 min(trimmean(DATA.RawData.SSCA,99)+std(DATA.RawData.SSCA)*2,2.5E5)];


axes(MyHandles.Ax(1));
     set(MyHandles.DAT(1).PTS,'XData',DATA.RawData.FSCH,'YData',DATA.RawData.SSCH);
     set(MyHandles.Ax(1),'xlim',FSCH_LIM,'ylim',SSCH_LIM);    
axes(MyHandles.Ax(2));
     set(MyHandles.DAT(2).PTS,'XData',DATA.RawData.FSCA,'YData',DATA.RawData.FSCH);
     set(MyHandles.Ax(2),'xlim',FSCA_LIM,'ylim',FSCH_LIM);
axes(MyHandles.Ax(3));
     set(MyHandles.DAT(3).PTS,'XData',DATA.RawData.SSCA,'YData',DATA.RawData.SSCH);
     set(MyHandles.Ax(3),'xlim',SSCA_LIM,'ylim',SSCH_LIM);    
end

function PlotF()
GUIDAT = get(gcf,'UserData');
DATA = GUIDAT.DATA;
MASK = GUIDAT.DATA.MASK>0;

figure;
TempMask = sum(MASK,2)>2;
subplot(2,1,1)
line('XData',DATA.RawData.YFP(TempMask),'YData',DATA.RawData.FRET(TempMask),'LineStyle','none','Marker','o','MarkerSize',0.5,'MarkerFaceColor',[0.2 0.2 0.2],'MarkerEdgeColor','none');
 set(gca,'XScale','log','YScale','log')
xlim([10 10^6])
ylim([1 10^4])
subplot(2,1,2)
line('XData',DATA.RawData.CFP(TempMask),'YData',DATA.RawData.FRET(TempMask),'LineStyle','none','Marker','o','MarkerSize',0.5,'MarkerFaceColor',[0.2 0.2 0.2],'MarkerEdgeColor','none');
set(gca,'XScale','log','YScale','log')

function Save
GUIDAT = get(gcf,'UserData');
DATA = GUIDAT.DATA;
MASK = GUIDAT.DATA.MASK>0;
TempMask = sum(MASK,2)>2;
[fname pathname] = uiputfile('*.ffa', 'Save Filtered Data',[GUIDAT.WorkingPath '\' GUIDAT.Filename(1:end-4)]);
fx = fieldnames(DATA.RawData);
for k = 1:length(fx)
    Temp = DATA.RawData.(fx{k});    
    DATA.FILTData.(fx{k}) = Temp(TempMask);
end
save([pathname '\' fname],'-struct','DATA','Header','FILTData','RawData');

function SaveROI
GUIDAT = get(gcf,'UserData');
DAT = GUIDAT.DATA;
[fname pathname] = uiputfile('*.ROI', 'Save ROI',GUIDAT.WorkingPath);
save([pathname '\' fname], '-struct','DAT','ROI')

function LoadROI
FigNum = gcf;
GUIDAT = get(FigNum,'UserData');
[fname pathname] = uigetfile('*.ROI', 'Load ROI',GUIDAT.WorkingPath);
ROI = load([pathname '\' fname],'-mat');
bReplace = 1;
if length(GUIDAT.DATA.ROI(1).Edge.X)>1
    RePlace = questdlg('Replace ROI?','Yes','No');
    switch RePlace
        case 'Yes'
            bReplace = 1;
        otherwise
            bReplace = 0;
    end
end

if bReplace
    GUIDAT.DATA.ROI = ROI.ROI;
end
set(FigNum,'UserData',GUIDAT)
FilterData()

set(GUIDAT.MyHandles.ChBR1,'value',1);

set(GUIDAT.MyHandles.ChBR2,'value',1);

set(GUIDAT.MyHandles.ChBR3,'value',1);
RePlot()


function AdvanceFile
FigNum = gcf;
GUIDAT = get(FigNum,'UserData');
MyFiles = dir([GUIDAT.WorkingPath '\*.fcs']);
MyFileNames = {MyFiles.name};
idx = find(strcmp(MyFileNames,GUIDAT.Filename));
if idx<length(MyFiles)
    FACSData = ReadFACSFile([GUIDAT.WorkingPath MyFileNames{idx+1}]);
    GUIDAT.Filename = MyFileNames{idx+1};
    GUIDAT.DATA.RawData= FACSData.Data;
    GUIDAT.DATA.Header = FACSData.Header;
    if isfield(GUIDAT.DATA,'MASK')
        GUIDAT.DATA = rmfield(GUIDAT.DATA,'MASK');
    end
    set(FigNum, 'UserData',GUIDAT);
    set(FigNum,'Name',[GUIDAT.WorkingPath MyFileNames{idx+1}]);
    FilterData();
    RePlot();
end


function RetreatFile
FigNum = gcf;
GUIDAT = get(FigNum,'UserData');
MyFiles = dir([GUIDAT.WorkingPath '\*.fcs']);
MyFileNames = {MyFiles.name};
idx = find(strcmp(MyFileNames,GUIDAT.Filename));
if idx>1
    FACSData = ReadFACSFile([GUIDAT.WorkingPath MyFileNames{idx-1}]);
    GUIDAT.Filename = MyFileNames{idx-1};
    GUIDAT.DATA.RawData= FACSData.Data;
    GUIDAT.DATA.Header = FACSData.Header;
    if isfield(GUIDAT.DATA,'MASK')
        GUIDAT.DATA = rmfield(GUIDAT.DATA,'MASK');
    end
    set(FigNum, 'UserData',GUIDAT);
    set(FigNum,'Name',[GUIDAT.WorkingPath MyFileNames{idx-1}]);
    FilterData();
    RePlot();
end

function keypresser(varargin)
keypressed = varargin{2};
switch keypressed.Key
    case 's'
        Save();
    case 'rightarrow'
        AdvanceFile();
    case 'leftarrow'
        RetreatFile();        
    case 'o'
        LoadROI();
    case 'd'
        SaveROI();
end

function ProcessAllFiles
    FigNum = gcf;
    GUIDAT = get(FigNum,'UserData');
    MyFiles = dir([GUIDAT.WorkingPath '\*.fcs']);
    MyFileNames = {MyFiles.name};
    ROI = GUIDAT.DATA.ROI;
    fprintf('There are a total of %d Files,\nProceeding to process all ... \nHang Tight!\n',length(MyFileNames))
    
    for idx = 1:length(MyFiles)
        CurFileName = MyFileNames{idx};
        fprintf('File %d / %d : %s\n--------------------------------------\n', idx, length(MyFileNames),CurFileName)
        TempDATA = ReadFACSFile([GUIDAT.WorkingPath '\' CurFileName]);
        MASK = zeros(size(TempDATA.Data.FSCH,1),3);
        if length(ROI(1).Edge.X)>2
            [IN ON] = inpolygon(TempDATA.Data.FSCH,TempDATA.Data.SSCH,ROI(1).Edge.X,ROI(1).Edge.Y);
            MASK(:,1) = (IN | ON);
        end
        if length(ROI(2).Edge.X)>2
            [IN ON] = inpolygon(TempDATA.Data.FSCA,TempDATA.Data.FSCH,ROI(2).Edge.X,ROI(2).Edge.Y);
            MASK(:,2) = (IN | ON);
        end
        if length(ROI(3).Edge.X)>2
            [IN ON] = inpolygon(TempDATA.Data.SSCA,TempDATA.Data.SSCH,ROI(3).Edge.X,ROI(3).Edge.Y);
            MASK(:,3) = (IN | ON);
        end     
        TempMask = sum(MASK(:,1:3),2)>2;
        fprintf('Filtering kept %5.2f%% data points including %d cells \n',sum(TempMask)/size(MASK,1)*100,sum(TempMask))
        fx = fieldnames(TempDATA.Data);
        for k = 1:length(fx)
            Temp = TempDATA.Data.(fx{k});    
            TempDATA.FILTData.(fx{k}) = Temp(TempMask);
        end
        
        save([GUIDAT.WorkingPath CurFileName(1:end-3) 'ffa'],'-struct','TempDATA','Header','FILTData');
        fprintf('Wrote Filtered File: %s \n\n',[CurFileName(1:end-3) 'ffa'])
    end
    
    