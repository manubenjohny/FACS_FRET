% FRET 2-hybrid assay implemented using a flow cytometer -- LSR II and BD FACS Diva
% Manu Ben-Johny, mbj2124@cumc.columbia.edu

function output = FACS_FRET(varargin)
if nargin==0
    InitializeGUI;
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        rethrow(lasterror);
    end
end

function InitializeGUI()
clear global DATASTRUCT
global DATASTRUCT
MyHandles.FG = figure('name',['Flow FRET Analysis'], 'NumberTitle','off','units','normalized', 'position', [0 0.3 0.8 0.6],'color', [1 1 1],'KeyPressFcn',@keypresser);
i = 1; W= 0.2; H = 0.15;S = 0.01;
MyHandles.Constructs        = CreateFileHandler([S S*i+(i-1)*H W H],'Constructs'); i = i+1;
MyHandles.Spurious_FRET     = CreateFileHandler([S S*i+(i-1)*H W H],'Spurious_FRET');i = i+1;
MyHandles.Dimers            = CreateFileHandler([S S*i+(i-1)*H W H],'Dimers');i = i+1;
MyHandles.Donor             = CreateFileHandler([S S*i+(i-1)*H W H],'Donor');i = i+1;
MyHandles.Acceptor          = CreateFileHandler([S S*i+(i-1)*H W H],'Acceptor');i = i+1;
MyHandles.Background        = CreateFileHandler([S S*i+(i-1)*H W H],'Background');i = i+1;
MyHandles.BasePath          = uicontrol('style','pushbutton', 'units','normalized','position',[S S*i+(i-1)*H W/4 0.02],'string','Pick Directory','CallBack',['FACS_FRET(''BasePath'')']);
MyHandles.BasePathText      = uicontrol('style','text', 'units','normalized','position',[S+W/4+S S*i+(i-1)*H 1-3*S-W/4 0.02],'string','','HorizontalAlignment','left','FontSize',9);
i = i-1;

PanelAxPos = [0.2200    0.4900    0.7700    0.2600];
MyHandles.PanelAx           = uipanel('Position',PanelAxPos,'BackgroundColor',[ 1 1 1]*0.97,'visible','on');
MyHandles.bBGSubtract       = uicontrol('parent',MyHandles.PanelAx,'style','Checkbox','BackgroundColor',[ 1 1 1]*0.97,'units','normalized','position',[0.85 .9 0.15 0.1],'CallBack',['FACS_FRET(''UpdatePlots'','''')'],'String','Subtract Background','FontSize',8);
MyHandles.bPlotAll          = uicontrol('parent',MyHandles.PanelAx,'style','Checkbox','BackgroundColor',[ 1 1 1]*0.97,'units','normalized','position',[0.65 .9 0.15 0.1],'CallBack',['FACS_FRET(''UpdatePlots'','''')'],'String','Plot All','FontSize',8,'value',0);
MyHandles.Ax(1) = axes('parent',MyHandles.PanelAx,'Units', 'normalized','OuterPosition', [0 0 .33333 0.95],'TickDir','out','box','off','XLimMode','manual'); xlabel('CFP'); ylabel('YFP')
MyHandles.Ax(2) = axes('parent',MyHandles.PanelAx,'Units', 'normalized','OuterPosition', [.33333 0 .33333 0.95],'TickDir','out','box','off','XLimMode','manual'); xlabel('CFP'); ylabel('FRET')
MyHandles.Ax(3) = axes('parent',MyHandles.PanelAx,'Units', 'normalized','OuterPosition', [.33333*2 0 .33333 0.95],'TickDir','out','box','off','XLimMode','manual'); xlabel('YFP'); ylabel('FRET')

PanelButtonPos = [0.2200    0.9200    0.7700    0.0300];
MyHandles.PanelButton = uipanel('Position',PanelButtonPos,'BackgroundColor',[ 1 1 1],'BorderType','none'); BGW = 1/14;BGS = 0*BGW/16; j = 1;
MyHandles.AnalysisPanel.BG = uicontrol('parent',MyHandles.PanelButton, 'style','togglebutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Background','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''SetupBackground'')'],'SelectionHighlight','off'); j = j+1;
MyHandles.AnalysisPanel.RA = uicontrol('parent',MyHandles.PanelButton, 'style','togglebutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Setup RA','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''SetupRA'')'],'SelectionHighlight','off');j =j+1;
MyHandles.AnalysisPanel.RD = uicontrol('parent',MyHandles.PanelButton, 'style','togglebutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Setup RD','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''SetupRD'')'],'SelectionHighlight','off');j =j+1;
MyHandles.AnalysisPanel.Dimers = uicontrol('parent',MyHandles.PanelButton, 'style','togglebutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Setup Dimers','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''SetupDimers'')'],'SelectionHighlight','off');j = j+1;
MyHandles.AnalysisPanel.SpuriousFRET = uicontrol('parent',MyHandles.PanelButton, 'style','togglebutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Spurious FRET','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''SetupSpurious'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.SelSubPop = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Filtering','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''Filtering'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.AutoFit = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','AutoFit','BackgroundColor',[1 0.7 0.7], 'CallBack',['FACS_FRET(''AutoFitCaller'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.SaveCalib = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Save Calibration','BackgroundColor',[0.8 0.8 1], 'CallBack',['FACS_FRET(''SaveCalibration'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.LoadCalib = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Load Calibration','BackgroundColor',[0.8 0.8 1], 'CallBack',['FACS_FRET(''LoadCalibration'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.SaveAnalysis = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Save Analysis','BackgroundColor',[0.8 0.8 1], 'CallBack',['FACS_FRET(''SaveAnalysis'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.LoadAnalysis = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Load Analysis','BackgroundColor',[0.8 0.8 1], 'CallBack',['FACS_FRET(''LoadAnalysis'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.ExportFigure = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Export Fig','BackgroundColor',[0.8 0.8 1], 'CallBack',['FACS_FRET(''ExportFigure'')'],'SelectionHighlight','off');j=j+1;
MyHandles.AnalysisPanel.ExportData = uicontrol('parent',MyHandles.PanelButton, 'style','pushbutton', 'units','normalized','position',[(BGW+BGS)*(j-1) 0 BGW 1],'string','Export Data','BackgroundColor',[0.8 0.8 1], 'CallBack',['FACS_FRET(''ExportData'')'],'SelectionHighlight','off');
MyHandles.AnalysisPanel.ShowHidePanel = uicontrol('parent',MyHandles.PanelButton, 'style','togglebutton', 'units','normalized','position',[1-BGW 0 BGW 1],'string','Raw Data','BackgroundColor',[0.7 1 0.7], 'CallBack',['FACS_FRET(''ShowHideRawData'')'],'SelectionHighlight','off','value',1);

i = i - 1;
PanelPARAMSColor = [ 0.98 0.98 1];
PanelPARAMS2Color = [ 0.95 1 0.95];
PanelParamsPos = [ 0.2200    0.770    0.55    0.1300];
PanelParams2Pos = [ 0.77    0.770    0.22    0.1300];
MyHandles.PanelPARAMS = uipanel('Position',PanelParamsPos,'BackgroundColor',PanelPARAMSColor);
MyHandles.PanelPARAMS2 = uipanel('Position',PanelParams2Pos,'BackgroundColor',PanelPARAMS2Color);

BGW = 1/11;BGS = BGW/15; j =1;
MyHandles.ParameterPanel.BGValsC.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','BG.CFP','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.BGValsY.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','BG.YFP','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.BGValsF.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.03 BGW 0.25],'string','BG.FRET','BackgroundColor',PanelPARAMSColor); j = j+1;
MyHandles.ParameterPanel.BGValsC.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',1)']);
MyHandles.ParameterPanel.BGValsY.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',2)']);
MyHandles.ParameterPanel.BGValsF.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.05 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',3)']);
j = j+1;
MyHandles.ParameterPanel.BGValsC_cut.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','BG.CFP_Cutoff','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.BGValsY_cut.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','BG.YFP_Cutoff','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.BGValsF_cut.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.03 BGW 0.25],'string','BG.FRET_Cutoff','BackgroundColor',PanelPARAMSColor); j = j+1;
MyHandles.ParameterPanel.BGValsC_cut.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',1)']);
MyHandles.ParameterPanel.BGValsY_cut.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',2)']);
MyHandles.ParameterPanel.BGValsF_cut.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.05 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',3)']);
j = j+1;
MyHandles.ParameterPanel.RA1.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','RA1','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.RD1.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','RD1','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.RD2.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.03 BGW 0.25],'string','RD2','BackgroundColor',PanelPARAMSColor); 
j=j+1;
MyHandles.ParameterPanel.RA1.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',4)']);
MyHandles.ParameterPanel.RD1.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',5)']);
MyHandles.ParameterPanel.RD2.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.05 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',6)']);
j = j+1;
MyHandles.ParameterPanel.G.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','Gratio','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.F.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','Fratio','BackgroundColor',PanelPARAMSColor);
j = j+1;
MyHandles.ParameterPanel.G.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',7)']);
MyHandles.ParameterPanel.F.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',8)']);
j = j+1;
MyHandles.ParameterPanel.EASlope.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','EAslope','BackgroundColor',PanelPARAMSColor);
MyHandles.ParameterPanel.EDSlope.Label = uicontrol('parent',MyHandles.PanelPARAMS, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','EDslope','BackgroundColor',PanelPARAMSColor);
j = j+1;
MyHandles.ParameterPanel.EASlope.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',9)']);
MyHandles.ParameterPanel.EDSlope.Edit = uicontrol('parent',MyHandles.PanelPARAMS, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',10)']);
BGW = 1/5;BGS = BGW/5; j =1;
MyHandles.ParameterPanel.Kdeff.Label = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','Kdeff','BackgroundColor',PanelPARAMS2Color);
MyHandles.ParameterPanel.nDnA.Label = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','nD/nA','BackgroundColor',PanelPARAMS2Color);
j = j+1;
MyHandles.ParameterPanel.Kdeff.Edit = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',11)']);
MyHandles.ParameterPanel.nDnA.Edit = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',14)']);
j = j+1;
MyHandles.ParameterPanel.EAmax.Label = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.67 BGW 0.25],'string','EAmax','BackgroundColor',PanelPARAMS2Color);
MyHandles.ParameterPanel.EDmax.Label = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','text', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.33 BGW 0.25],'string','EDmax','BackgroundColor',PanelPARAMS2Color);
j=j+1;
MyHandles.ParameterPanel.EAmax.Edit = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.7 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',12)']);
MyHandles.ParameterPanel.EDmax.Edit = uicontrol('parent',MyHandles.PanelPARAMS2, 'style','edit', 'units','normalized','position',[(BGW+BGS)*(j-1) 0.35 BGW 0.25],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateCalib'',13)']);

DATASTRUCT.ActivePath = '';

i = i-3;
PanelAx2Pos = [0.2200    0.0100    0.7700    0.4550];
MyHandles.PanelAx2 = uipanel('Position',[2*S+W S*(i-1)+(i-2)*H 1-3*S-W 3.3*H-4*S],'BackgroundColor',[ 1 1 1]);
MyHandles.Ax(4) = axes('parent',MyHandles.PanelAx2,'Units', 'normalized','OuterPosition', [0 0 .5 0.95],'TickDir','out','box','off','XLimMode','manual'); xlabel('Dfree'); ylabel('E_A')
MyHandles.Ax(5) = axes('parent',MyHandles.PanelAx2,'Units', 'normalized','OuterPosition', [0.5 0 .5 0.95],'TickDir','out','box','off','XLimMode','manual'); xlabel('Afree'); ylabel('E_D')

% Setup menus
MyHandles.hcmenu = zeros(1,length(MyHandles.Ax));
for kj = 1:length(MyHandles.Ax)
    MyHandles.hcmenu(kj) = uicontextmenu;
    MyHandles.CM(kj) = uimenu(MyHandles.hcmenu(kj), 'Label', 'Set Limits', 'Callback', ['FACS_FRET(''ChangeAxLimits'',' num2str(kj) ')']);
    set(MyHandles.Ax(kj),'uicontextmenu',MyHandles.hcmenu(kj));
end


%initialize Parameters
DATASTRUCT.PARAMETERS.BG.CFP = 0;
DATASTRUCT.PARAMETERS.BG.YFP = 0;
DATASTRUCT.PARAMETERS.BG.FRET = 0;
DATASTRUCT.PARAMETERS.RA1 = 0;
DATASTRUCT.PARAMETERS.RD1 = 0;
DATASTRUCT.PARAMETERS.RD2 = 0;
DATASTRUCT.PARAMETERS.gRatio = 0;
DATASTRUCT.PARAMETERS.fRatio = 0;

DATASTRUCT.PARAMETERS.EASlope = 0;
DATASTRUCT.PARAMETERS.EDSlope = 0;

DATASTRUCT.PARAMETERS.BG.CFP_CUTOFF  = 0;
DATASTRUCT.PARAMETERS.BG.YFP_CUTOFF = 0;
DATASTRUCT.PARAMETERS.BG.FRET_CUTOFF = 0;

DATASTRUCT.PARAMETERS.YFP_MAX = 4E8;
DATASTRUCT.PARAMETERS.CFP_MAX = 4E8;
DATASTRUCT.PARAMETERS.FRET_MAX = 4E8;

DATASTRUCT.WORKFLOW.bRASet = 0;
DATASTRUCT.WORKFLOW.bRDSet = 0;
DATASTRUCT.WORKFLOW.bBGSet = 0;
DATASTRUCT.WORKFLOW.bDimerSet = 0;

DATASTRUCT.WORKFLOW.bRASet = 0;
DATASTRUCT.WORKFLOW.bRDSet = 0;
DATASTRUCT.WORKFLOW.bBGSet = 0;
DATASTRUCT.WORKFLOW.bDimerSet = 0;
DATASTRUCT.WORKFLOW.bMaskSet = 0;


DATASTRUCT.MASKPARAMS.CFP.Min = DATASTRUCT.PARAMETERS.BG.CFP_CUTOFF;
DATASTRUCT.MASKPARAMS.CFP.Max = 2.5E5;
DATASTRUCT.MASKPARAMS.YFP.Min = DATASTRUCT.PARAMETERS.BG.YFP_CUTOFF;
DATASTRUCT.MASKPARAMS.YFP.Max = 2.5E5;
DATASTRUCT.MASKPARAMS.FRET.Min = DATASTRUCT.PARAMETERS.BG.FRET_CUTOFF;
DATASTRUCT.MASKPARAMS.FRET.Max = 2.5E5;
DATASTRUCT.MASKPARAMS.NDNA = 50;
DATASTRUCT.MASKPARAMS.NAND = 50;

setPARAMETERS(DATASTRUCT.PARAMETERS, MyHandles)

set(MyHandles.BasePathText,'string',DATASTRUCT.ActivePath);

DATASTRUCT.Handles = MyHandles;
% set(MyHandles.FG,'userdata',DATASTRUCT)

function MyHandles = CreateFileHandler(Position,Tag)
MyHandles.MainPanel = uipanel('Title',Tag,'position',Position,'BackgroundColor',[1 1 1]*0.95,'BorderType','etchedout','ForegroundColor',[1 1 1]*0,'FontSize',10);
MyHandles.Plus = uicontrol('style','pushbutton','parent',MyHandles.MainPanel,'units','normalized','Position',[0.9 0.9 0.1 0.2],'String','+','BackgroundColor',[1 1 1]*0,'ForegroundColor',[1 1 1],'FontSize',12,'FontWeight','bold','CallBack',['FACS_FRET(''AddFile'',''' Tag ''')']);
MyHandles.Minus = uicontrol('style','pushbutton','parent',MyHandles.MainPanel,'units','normalized','Position',[0.8 .9 0.1 0.2],'String','-','BackgroundColor',[1 1 1]*0,'ForegroundColor',[1 1 1],'FontSize',12,'FontWeight','bold','CallBack',['FACS_FRET(''RemoveFile'',''' Tag ''')']);
MyHandles.FileList = uicontrol('style','listbox','parent',MyHandles.MainPanel,'units','normalized','Position',[0.01 .01 0.98 0.85], 'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdatePlots'',''' Tag ''')'],'string',{});

function Handles = AddFile(Tag)
%     DATASTRUCT = get(gcf,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    if ~isempty(DATASTRUCT.ActivePath)
        [FNAME PATH] = uigetfile('*.ffa', ['Choose ' Tag ' File'],DATASTRUCT.ActivePath);
    else
        [FNAME PATH] = uigetfile('*.ffa', ['Choose ' Tag ' File']);
    end
    FileList = get(MyHandles.(Tag).FileList,'string');
    nFiles = length(FileList);
    TempFile = [PATH FNAME];    
    if length(TempFile)<3
        return;
    end
    if isempty(FileList) || (sum(strcmp(FileList,TempFile))<1)
        FileList{nFiles+1} = TempFile;
        set(MyHandles.(Tag).FileList,'string',FileList,'value',nFiles+1);
    end
    DATASTRUCT.(Tag).File(nFiles+1) = loadFFA(TempFile); 
    UpdatePlots(Tag)
   
function Handles = BasePath()
global DATASTRUCT;
% DATASTRUCT = get(gcf,'userdata');
MyHandles = DATASTRUCT.Handles;
DATASTRUCT.ActivePath = uigetdir;
if DATASTRUCT.ActivePath ==0
    DATASTRUCT.ActivePath = '';
end
set(MyHandles.BasePathText,'string',DATASTRUCT.ActivePath);
% set(gcf,'userdata',DATASTRUCT);

function Handles = RemoveFile(Tag)
%     DATASTRUCT = get(gcf,'userdata');
    global DATASTRUCT;
    MyHandles = DATASTRUCT.Handles;
    FileList = get(MyHandles.(Tag).FileList,'string');
    ListX = 1:length(FileList);
    sel = get(MyHandles.(Tag).FileList,'value');
    File2Delete = FileList{sel};
    bDelete = questdlg(['Remove File ' FileList{sel}],'Flow FRET Analysis');    
    switch bDelete
        case 'Yes'
            FileList = FileList(~(ListX==sel));
            set(MyHandles.(Tag).FileList,'string',FileList,'value',max(sel-1,1));
            DATASTRUCT.(Tag).File = DATASTRUCT.(Tag).File(~(ListX==sel));   
            if isempty(DATASTRUCT.(Tag).File)
                DATASTRUCT = rmfield(DATASTRUCT,Tag);
            end
    end
%     set(gcf,'userdata',DATASTRUCT);
    UpdatePlots(Tag)

function FCSData = loadFFA(MyFile)
    Temp = load(MyFile,'-mat','Header','FILTData'); 
    FCSData.Data = Temp.FILTData;
    FitParams.Kdeff = 0;
    FitParams.EAmax = 0;
    FitParams.EDmax = 0;
    FitParams.nDnA  = 1;
    FCSData.FitParams = FitParams;
    
    
    
function UpdatePlots(Tag)
    global DATASTRUCT;
%     DATASTRUCT = get(gcf,'userdata');
    MyHandles = DATASTRUCT.Handles;
    PARAMS = getPARAMETERS(MyHandles);

    if isempty(Tag)
        Tag = get(MyHandles.Ax(1),'UserData'); 
    end
    if ~isfield(DATASTRUCT, Tag)
        ClearAllAxes();
        return;
    end
    if ~isempty(Tag)
         if get(MyHandles.bPlotAll,'value') < 1
             DATA2Plot = DATASTRUCT.(Tag).File(get(MyHandles.(Tag).FileList,'value'));
         else
             DATA2Plot = DATASTRUCT.(Tag).File;
         end
        switch Tag
            case 'Background'
                bSubtract = 0;
            otherwise;
                bSubtract = get(MyHandles.bBGSubtract,'value');
            end
                 
         if length(DATA2Plot)>1
            MyColors = jet(length(DATA2Plot));
         else
             MyColors = [0 0 0];
         end 
         cla(MyHandles.Ax(1))               
         cla(MyHandles.Ax(2))
         cla(MyHandles.Ax(3))
         
         CFP_LIMIT = [0 0];
         YFP_LIMIT = [0 0];
         FRET_LIMIT = [0 0];
         for i = 1:length(DATA2Plot)
             
             MyCFP = DATA2Plot(i).Data.CFP - bSubtract.*PARAMS.BG.CFP;
             MyYFP = DATA2Plot(i).Data.YFP - bSubtract.*PARAMS.BG.YFP;
             MyFRET = DATA2Plot(i).Data.FRET - bSubtract.*PARAMS.BG.FRET;
             
             MASK = ones(size(MyCFP))>0;
             
             CFP_LIMIT(1) = min([min(CFP_LIMIT) min(MyCFP)]);
             CFP_LIMIT(2) = max([max(CFP_LIMIT) prctile(MyCFP,99.9)]);
             YFP_LIMIT(1) = min([min(YFP_LIMIT) min(MyYFP)]);
             YFP_LIMIT(2) = max([max(YFP_LIMIT) prctile(MyYFP,99.9)]);
             FRET_LIMIT(1) = min([min(FRET_LIMIT) min(MyFRET)]);
             FRET_LIMIT(2) = max([max(FRET_LIMIT) prctile(MyFRET,99.9)]);
             if get(MyHandles.AnalysisPanel.ShowHidePanel,'value')>0
                 axes(MyHandles.Ax(1))
                 line(MyCFP(MASK),MyYFP(MASK),'Marker','o','LineStyle','none','MarkerSize',1,'MarkerEdgeColor','none','MarkerFaceColor',MyColors(i,:)) 
                 set(MyHandles.Ax(1),'UserData',Tag)    
                 if strcmp(get(MyHandles.CM(1),'checked'),'off')
                    xlim(CFP_LIMIT); ylim(YFP_LIMIT);
                 end

                 axes(MyHandles.Ax(2))
                 line(MyCFP(MASK),MyFRET(MASK),'Marker','o','LineStyle','none','MarkerSize',1,'MarkerEdgeColor','none','MarkerFaceColor',MyColors(i,:))
                 set(MyHandles.Ax(2),'UserData',Tag)
                 if strcmp(get(MyHandles.CM(2),'checked'),'off')
                    xlim(CFP_LIMIT); ylim(FRET_LIMIT);
                 end
                 
                 axes(MyHandles.Ax(3))
                 line(MyYFP(MASK),MyFRET(MASK),'Marker','o','LineStyle','none','MarkerSize',1,'MarkerEdgeColor','none','MarkerFaceColor',MyColors(i,:))
                 set(MyHandles.Ax(3),'UserData',Tag)
                 if strcmp(get(MyHandles.CM(3),'checked'),'off')
                    xlim(YFP_LIMIT); ylim(FRET_LIMIT);
                 end
             end
         end
         switch Tag
             case 'Background'
                 
                 if ~isempty(PARAMS.BG.CFP) && ~isempty(PARAMS.BG.YFP) && ~isempty(PARAMS.BG.FRET)
                     axes(MyHandles.Ax(1))
                     line(PARAMS.BG.CFP,PARAMS.BG.YFP,'Marker','o','LineStyle','none','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]) 
                     set(MyHandles.Ax(1),'UserData',Tag)    
                     if strcmp(get(MyHandles.CM(1),'checked'),'off')
                        xlim(CFP_LIMIT); ylim(YFP_LIMIT);
                     end
                     
                     axes(MyHandles.Ax(2))
                     line(PARAMS.BG.CFP,PARAMS.BG.FRET,'Marker','o','LineStyle','none','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]) 
                     set(MyHandles.Ax(2),'UserData',Tag)
                     if strcmp(get(MyHandles.CM(2),'checked'),'off')
                        xlim(CFP_LIMIT); ylim(FRET_LIMIT);
                     end
                     axes(MyHandles.Ax(3))
                     line(PARAMS.BG.YFP,PARAMS.BG.FRET,'Marker','o','LineStyle','none','MarkerSize',8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]) 
                     set(MyHandles.Ax(3),'UserData',Tag)
                     if strcmp(get(MyHandles.CM(3),'checked'),'off')
                        xlim(YFP_LIMIT); ylim(FRET_LIMIT);
                     end
                 end
            case 'Acceptor' 

                 RA = str2num(get(MyHandles.ParameterPanel.RA1.Edit,'string'));
                 if ~isempty(RA)
                     FRET_FIT = RA*[0 max(YFP_LIMIT)];
                     axes(MyHandles.Ax(3))
                     line([0 max(YFP_LIMIT)], FRET_FIT, 'marker','none','LineStyle','-', 'LineWidth',2,'Color',[0 1 0]);
                     if strcmp(get(MyHandles.CM(3),'checked'),'off')
                        xlim(YFP_LIMIT); ylim(FRET_LIMIT);
                     end
                 end
                 
            case 'Donor' 
                 set(MyHandles.bBGSubtract,'value',1)
                 RD1 = str2num(get(MyHandles.ParameterPanel.RD1.Edit,'string'));                 
                 RD2 = str2num(get(MyHandles.ParameterPanel.RD2.Edit,'string'));
                 if ~isempty(RD1) && ~isempty(RD2)
                     FRET_FIT = RD1*[0 max(CFP_LIMIT)];
                     YFP_FIT = RD2*[0 max(CFP_LIMIT)];
                     axes(MyHandles.Ax(2))
                        line([0 max(CFP_LIMIT)],FRET_FIT,'Marker','none','LineStyle','-','Color',[0 0 1],'LineWidth',2);
                        if strcmp(get(MyHandles.CM(2),'checked'),'off')
                            xlim(CFP_LIMIT); ylim(FRET_LIMIT);
                        end
                     axes(MyHandles.Ax(1))
                        line([0 max(CFP_LIMIT)],YFP_FIT,'Marker','none','LineStyle','-','Color',[0 0 1],'LineWidth',2);
                        if strcmp(get(MyHandles.CM(1),'checked'),'off')
                            xlim(CFP_LIMIT); ylim(YFP_LIMIT);   
                        end
                 end
             case 'Dimers'
                    if get(MyHandles.bPlotAll,'value')==1 
                        sel =1:length(get(MyHandles.(Tag).FileList,'string'));
                    else
                        sel = get(MyHandles.(Tag).FileList,'value');
                    end
                    SUBDATA = DATASTRUCT.('Dimers').File;
                    cla(MyHandles.Ax(4));
                    cla(MyHandles.Ax(5));   
                    Gx = str2num(get(MyHandles.ParameterPanel.G.Edit,'string'));
                    Fx = str2num(get(MyHandles.ParameterPanel.F.Edit,'string'));
                    DATASTRUCT.PARAMETERS.fRatio = Fx;
                    DATASTRUCT.PARAMETERS.gRatio = Gx;
                    DATASTRUCT.PARAMETERS.BRatio = Fx*Gx;
                    
                    
                    if DATASTRUCT.WORKFLOW.bDimerSet>0
                        EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
                        EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
                        if isempty(EASlope) 
                            EASlope = 0;
                        end
                        if isempty(EDSlope) 
                            EDSlope = 0;
                        end
                        
                        NDNAmean = zeros(1,length(SUBDATA));
                        
                        for i = sel
                          DATASTRUCT.('Dimers').File(i).Data.EA = DATASTRUCT.('Dimers').File(i).Data.YFPfret./DATASTRUCT.('Dimers').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
                          DATASTRUCT.('Dimers').File(i).Data.ED = DATASTRUCT.('Dimers').File(i).Data.YFPfret./(DATASTRUCT.('Dimers').File(i).Data.YFPfret+ DATASTRUCT.('Dimers').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);

                          ND = DATASTRUCT.('Dimers').File(i).Data.CFPdirect./(1-DATASTRUCT.('Dimers').File(i).Data.ED);
                          NA = DATASTRUCT.('Dimers').File(i).Data.YFPdirect/(DATASTRUCT.PARAMETERS.fRatio*DATASTRUCT.PARAMETERS.gRatio);
                          
                          NDNAtemp  =ND./NA; 
                          NDNAmean(i) = median(NDNAtemp(~isnan(NDNAtemp)));
                          

                          DATASTRUCT.('Dimers').File(i).Data.EA = DATASTRUCT.('Dimers').File(i).Data.EA - ND*EASlope;
                          DATASTRUCT.('Dimers').File(i).Data.ED = DATASTRUCT.('Dimers').File(i).Data.ED - NA*EDSlope;

                          Cutoff = 1;
                          CFPTEMP  = [DATASTRUCT.('Dimers').File(i).Data.CFP]  - PARAMS.BG.CFP;
                          YFPTEMP  = [DATASTRUCT.('Dimers').File(i).Data.YFP]  - PARAMS.BG.YFP;
                          FRETTEMP = [DATASTRUCT.('Dimers').File(i).Data.FRET] - PARAMS.BG.FRET;
                          
                          MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
                          MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
            
                          if length(sel)>1
                              TempColor = MyColors(i,:);
                          else
                              TempColor = MyColors;
                          end
                          axes(MyHandles.Ax(4));
                          line(ND(MASK) ,DATASTRUCT.('Dimers').File(i).Data.EA(MASK),'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor',TempColor);
                          if strcmp(get(MyHandles.CM(4),'checked'),'off')
                            ylim([-0.1 0.8]);xlim([0 4E4])
                          end
                          axes(MyHandles.Ax(5)); 
                          line(NA(MASK),DATASTRUCT.('Dimers').File(i).Data.ED(MASK),'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor',TempColor);
                          if strcmp(get(MyHandles.CM(5),'checked'),'off')
                            ylim([-0.1 0.8]);xlim([0 4E4])
                          end
                        end
                        fprintf('Dimer Summary \n---------------------------------------------------\n')
                        TempFiles = get(MyHandles.(Tag).FileList,'string');
                        for kj = 1:length(sel)
                         fprintf('%s \t\t %5.4f\n',TempFiles{sel(kj)},NDNAmean(sel(kj)))
                        end
                        fprintf('---------------------------------------------------\n')
%                           figure; 
%                           bar(NDNAmean)
%                           ylim([0 2])



                    end
             case 'Spurious_FRET'
                 SetupSpurious()
                 
                 
             case 'Constructs'
                 SetupConstructs();
                    
         end
    end
    set(MyHandles.Ax(1),'UserData',Tag)
    
function DAT = mergeData(SubData)   
    numFiles = length(SubData);
    nPts = 0;
    for i = 1:numFiles
        nPts = nPts+length(SubData(i).Data.CFP);
    end
    init = 0;
    DAT.CFP = zeros(1,nPts);
    DAT.YFP = DAT.CFP;
    DAT.FRET = DAT.CFP;    
    for i = 1:numFiles
        DAT.CFP(init+1:init+length(SubData(i).Data.CFP)) = SubData(i).Data.CFP;
        DAT.YFP(init+1:init+length(SubData(i).Data.CFP)) = SubData(i).Data.YFP;
        DAT.FRET(init+1:init+length(SubData(i).Data.CFP)) = SubData(i).Data.FRET;
        init = init+length(SubData(i).Data.CFP);
    end
    
    DAT.CFP = DAT.CFP';
    DAT.YFP = DAT.YFP';
    DAT.FRET = DAT.FRET';
    
function SetupBackground(varargin)
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    if nargin == 0
        bSuppress = 0; 
    else
        bSuppress = varargin{1}>0;
    end
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    if  (DATASTRUCT.WORKFLOW.bBGSet == 0);   
        % Check if there is any data for Background. 
        if ~isfield(DATASTRUCT,'Background')
            set(MyHandles.AnalysisPanel.BG,'Value',0);
            return;
        end
        SUBDATA = DATASTRUCT.('Background').File;
        if ~isempty(SUBDATA) 
            DAT = mergeData(SUBDATA);
            CFPTEMP  = DAT.CFP;
            YFPTEMP  = DAT.YFP;
            FRETTEMP = DAT.FRET;
                        
            DATASTRUCT.PARAMETERS.BG.CFP  = trimmean(CFPTEMP,0.1);
            DATASTRUCT.PARAMETERS.BG.YFP = trimmean(YFPTEMP,0.1);
            DATASTRUCT.PARAMETERS.BG.FRET = trimmean(FRETTEMP,0.1);
            DATASTRUCT.PARAMETERS.BG.CFP_CUTOFF  = prctile(CFPTEMP,99.5) - DATASTRUCT.PARAMETERS.BG.CFP;
            DATASTRUCT.PARAMETERS.BG.YFP_CUTOFF = prctile(YFPTEMP,99.5)  - DATASTRUCT.PARAMETERS.BG.YFP;
            DATASTRUCT.PARAMETERS.BG.FRET_CUTOFF = prctile(FRETTEMP,99.5)- DATASTRUCT.PARAMETERS.BG.FRET;
            
            DATASTRUCT.MASKPARAMS.CFP.Min = DATASTRUCT.PARAMETERS.BG.CFP_CUTOFF;
            DATASTRUCT.MASKPARAMS.YFP.Min = DATASTRUCT.PARAMETERS.BG.YFP_CUTOFF;
            DATASTRUCT.MASKPARAMS.FRET.Min = DATASTRUCT.PARAMETERS.BG.FRET_CUTOFF;
            
        end
        set(MyHandles.AnalysisPanel.BG,'BackgroundColor',[0.7 1 0.7],'Value',1);
        DATASTRUCT.WORKFLOW.bBGSet = 1;    
    else 
        set(MyHandles.AnalysisPanel.BG,'BackgroundColor',[1 0.7 0.7])
        set(MyHandles.AnalysisPanel.BG,'Value',0);
        
        
        DATASTRUCT.PARAMETERS.BG.CFP = 0;
        DATASTRUCT.PARAMETERS.BG.YFP = 0;
        DATASTRUCT.PARAMETERS.BG.FRET = 0;   
        DATASTRUCT.PARAMETERS.BG.CFP_CUTOFF  = 0;
        DATASTRUCT.PARAMETERS.BG.YFP_CUTOFF = 0;
        DATASTRUCT.PARAMETERS.BG.FRET_CUTOFF = 0;
        
        DATASTRUCT.WORKFLOW.bBGSet = 0;
    end
    set(MyHandles.ParameterPanel.BGValsC.Edit,'string',DATASTRUCT.PARAMETERS.BG.CFP)
    set(MyHandles.ParameterPanel.BGValsY.Edit,'string',DATASTRUCT.PARAMETERS.BG.YFP)
    set(MyHandles.ParameterPanel.BGValsF.Edit,'string',DATASTRUCT.PARAMETERS.BG.FRET)
    setPARAMETERS(DATASTRUCT.PARAMETERS, MyHandles);
    if ~bSuppress
        UpdatePlots('Background')
    end
%     set(MainFig,'userdata',DATASTRUCT)
    
function SetupRA(varargin)
    if nargin == 0
        bSuppress = 0; 
    else
        bSuppress = varargin{1}>0;
    end
    global DATASTRUCT
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    MyHandles = DATASTRUCT.Handles;
    SUBDATA = DATASTRUCT.('Acceptor').File;
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    
    
    
    if ~isempty(SUBDATA) && (DATASTRUCT.WORKFLOW.bRASet==0);
        DAT = mergeData(SUBDATA);
        YRaw = DAT.YFP;
        FRaw = DAT.FRET;
        YFPTEMP  = YRaw - BG.YFP;
        FRETTEMP = FRaw - BG.FRET;
        MASK     = (YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
        YFPTEMP  = YFPTEMP(MASK);
        FRETTEMP = FRETTEMP(MASK);
        
        
        RA1 = iterativeScalarFit(YFPTEMP,FRETTEMP);
        
        DATASTRUCT.WORKFLOW.bRASet = 1;
        set(MyHandles.AnalysisPanel.RA,'BackgroundColor',[0.7 1 0.7],'Value',1);
        DATASTRUCT.PARAMETERS.RA1 = RA1;
        
        
        
    else
        DATASTRUCT.PARAMETERS.RA1 = 0;
        set(MyHandles.AnalysisPanel.RA,'BackgroundColor',[1 0.7 0.7])
        set(MyHandles.AnalysisPanel.RA,'Value',0);
        DATASTRUCT.WORKFLOW.bRASet = 0;
        RA1 = 0;
    end
    set(MyHandles.ParameterPanel.RA1.Edit,'string',RA1)
    if ~bSuppress
        UpdatePlots('Acceptor');
    end
%     set(MainFig,'userdata',DATASTRUCT)
   
function SetupRD(varargin)
    if nargin == 0
        bSuppress = 0; 
    else
        bSuppress = varargin{1}>0;
    end
    global DATASTRUCT
% %     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    MyHandles = DATASTRUCT.Handles;
    SUBDATA = DATASTRUCT.('Donor').File;     
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    if ~isempty(SUBDATA) && (DATASTRUCT.WORKFLOW.bRDSet == 0);
        DAT = mergeData(SUBDATA);
        YRaw = DAT.YFP;
        FRaw = DAT.FRET;
        CRaw = DAT.CFP;
        
        CFPTEMP  = [CRaw] - BG.CFP;
        YFPTEMP  = [YRaw] - BG.YFP;
        FRETTEMP = [FRaw] - BG.FRET;
        
        MASK = (CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
        CFPTEMP = CFPTEMP(MASK);
        YFPTEMP = YFPTEMP(MASK);
        FRETTEMP = FRETTEMP(MASK);
        RD1 = iterativeScalarFit(CFPTEMP,FRETTEMP);             
        RD2 = iterativeScalarFit(CFPTEMP, YFPTEMP);
        
        set(MyHandles.AnalysisPanel.RD,'BackgroundColor',[0.7 1 0.7],'Value',1);
        DATASTRUCT.WORKFLOW.bRDSet = 1;        
        DATASTRUCT.PARAMETERS.RD1 = RD1;
        DATASTRUCT.PARAMETERS.RD2 = RD2;
    else
        DATASTRUCT.PARAMETERS.RD1 = 0; RD1 = 0;       
        DATASTRUCT.PARAMETERS.RD2 = 0; RD2 = 0;
        DATASTRUCT.WORKFLOW.bRDSet = 0;
        set(MyHandles.AnalysisPanel.RD,'BackgroundColor',[1 0.7 0.7])
        set(MyHandles.AnalysisPanel.RD,'Value',0);
    end
    set(MyHandles.ParameterPanel.RD1.Edit,'string',RD1)
    set(MyHandles.ParameterPanel.RD2.Edit,'string',RD2)
    if ~bSuppress
    UpdatePlots('Donor');
    end
%     set(MainFig,'userdata',DATASTRUCT)
    
function SetupDimers(varargin)
    if nargin == 0
        bSuppress = 0; 
    else
        bSuppress = varargin{1}>0;
    end
     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    SUBDATA = DATASTRUCT.('Dimers').File;
    MyColors = jet(length(SUBDATA));
    
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;

    for i = 1:length(SUBDATA)
        CFPTEMP  = [DATASTRUCT.('Dimers').File(i).Data.CFP]  - BG.CFP;
        YFPTEMP  = [DATASTRUCT.('Dimers').File(i).Data.YFP]  - BG.YFP;
        FRETTEMP = [DATASTRUCT.('Dimers').File(i).Data.FRET] - BG.FRET;
        
        Cutoff = 1;
        MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
        MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
                 
        Fin = [CFPTEMP YFPTEMP FRETTEMP]';
        Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';
        DATASTRUCT.('Dimers').File(i).Data.CFPdirect  = Sout(:,1);
        DATASTRUCT.('Dimers').File(i).Data.YFPdirect  = Sout(:,2);
        DATASTRUCT.('Dimers').File(i).Data.YFPfret = Sout(:,3);
        
        DATASTRUCT.('Dimers').File(i).Data.CFPdirect(MASK<1) = 0;
        DATASTRUCT.('Dimers').File(i).Data.YFPdirect(MASK<1) = 0;
        DATASTRUCT.('Dimers').File(i).Data.YFPfret(MASK<1) = 0;
    end
    if ~bSuppress
        FigNew = figure('color',[ 1 1 1], 'name', 'Estimating G / F factors ...', 'units', 'normalized'); 
        getFRETStd(DATASTRUCT.('Dimers'), MainFig,FigNew);
        waitfor(FigNew)
    else
        
        SUBDATA = DATASTRUCT.('Dimers');
        CALIBDATA_median = zeros(length(SUBDATA.File),2);
        for kk = 1:length(SUBDATA.File)
            Rx = SUBDATA.File(kk).Data.YFPdirect./SUBDATA.File(kk).Data.CFPdirect;
            Ry = SUBDATA.File(kk).Data.YFPfret./SUBDATA.File(kk).Data.CFPdirect;
            MASKx = SUBDATA.File(kk).Data.CFPdirect==0;
            Rx = Rx(~MASKx);
            Ry = Ry(~MASKx);
            [Rx Ry] = removeNaN(Rx,Ry);            
            CALIBDATA_median(kk,:) = [median(Rx), median(Ry)];
        end
        b = robustfit(CALIBDATA_median(:,1),CALIBDATA_median(:,2),'fair',1.4);      
        set(MyHandles.ParameterPanel.G.Edit,'string',num2str(1/b(2)))
        set(MyHandles.ParameterPanel.F.Edit,'string',num2str(-b(1))) 
    end
    
    Gx = str2num(get(MyHandles.ParameterPanel.G.Edit,'string'));
    Fx = str2num(get(MyHandles.ParameterPanel.F.Edit,'string'));
    DATASTRUCT.PARAMETERS.fRatio = Fx;
    DATASTRUCT.PARAMETERS.gRatio = Gx;
    DATASTRUCT.PARAMETERS.BRatio = DATASTRUCT.PARAMETERS.fRatio*DATASTRUCT.PARAMETERS.gRatio;
    set(MyHandles.AnalysisPanel.Dimers,'BackgroundColor',[0.7 1 0.7]);
    DATASTRUCT.WORKFLOW.bDimerSet = 1;
   
%     set(MainFig,'userdata',DATASTRUCT)
    if ~bSuppress
        UpdatePlots('Dimers');
    else
        SUBDATA = DATASTRUCT.('Dimers');
        for ik = 1:length(SUBDATA.File)
          DATASTRUCT.('Dimers').File(ik).Data.EA = DATASTRUCT.('Dimers').File(ik).Data.YFPfret./DATASTRUCT.('Dimers').File(ik).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
          DATASTRUCT.('Dimers').File(ik).Data.ED = DATASTRUCT.('Dimers').File(ik).Data.YFPfret./(DATASTRUCT.('Dimers').File(ik).Data.YFPfret+ DATASTRUCT.('Dimers').File(ik).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);

          ND = DATASTRUCT.('Dimers').File(ik).Data.CFPdirect./(1-DATASTRUCT.('Dimers').File(ik).Data.ED);
          NA = DATASTRUCT.('Dimers').File(ik).Data.YFPdirect/(DATASTRUCT.PARAMETERS.fRatio*DATASTRUCT.PARAMETERS.gRatio);

            Cutoff = 1;
            CFPTEMP  = [DATASTRUCT.('Dimers').File(ik).Data.CFP]  - DATASTRUCT.PARAMETERS.BG.CFP;
            YFPTEMP  = [DATASTRUCT.('Dimers').File(ik).Data.YFP]  - DATASTRUCT.PARAMETERS.BG.YFP;
            FRETTEMP = [DATASTRUCT.('Dimers').File(ik).Data.FRET] - DATASTRUCT.PARAMETERS.BG.FRET;
          
          MASKt = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
          MASKt = MASKt&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
          
          NDNAtemp  =ND./NA; 
          NDNAtemp = NDNAtemp(MASKt);
          NDNAmean(ik) = median(NDNAtemp(~isnan(NDNAtemp)));      
        end
        fprintf('Dimer Summary \n---------------------------------------------------\n')
        TempFiles = get(MyHandles.('Dimers').FileList,'string');
        for kj = 1:length(TempFiles)
         fprintf('%s \t\t %5.4f\n',TempFiles{kj},NDNAmean(kj))
        end
        fprintf('---------------------------------------------------\n')
    end
   

function SetupSpurious()
    global DATASTRUCT
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    MyHandles = DATASTRUCT.Handles;
    
    SUBDATA = DATASTRUCT.('Spurious_FRET').File;
    sel = get(MyHandles.Spurious_FRET.FileList,'value');
    
    
    MyColors = 0*jet(length(SUBDATA));

    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;
    

    
    for i = sel

        CFPTEMP  = [DATASTRUCT.('Spurious_FRET').File(i).Data.CFP] - DATASTRUCT.PARAMETERS.BG.CFP;
        YFPTEMP  = [DATASTRUCT.('Spurious_FRET').File(i).Data.YFP] - DATASTRUCT.PARAMETERS.BG.YFP;
        FRETTEMP = [DATASTRUCT.('Spurious_FRET').File(i).Data.FRET] - DATASTRUCT.PARAMETERS.BG.FRET;
        
        MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
        MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);

        Fin = [CFPTEMP YFPTEMP FRETTEMP]';
        Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';
        
        DATASTRUCT.('Spurious_FRET').File(i).Data.CFPdirect  = Sout(:,1);
        DATASTRUCT.('Spurious_FRET').File(i).Data.YFPdirect  = Sout(:,2);
        DATASTRUCT.('Spurious_FRET').File(i).Data.YFPfret = Sout(:,3);

        DATASTRUCT.('Spurious_FRET').File(i).Data.CFPdirect(MASK<1) = 0;
        DATASTRUCT.('Spurious_FRET').File(i).Data.YFPdirect(MASK<1) = 0;
        DATASTRUCT.('Spurious_FRET').File(i).Data.YFPfret(MASK<1) = 0;

    end
    

    cla(MyHandles.Ax(4));
    cla(MyHandles.Ax(5));   
    for i = sel
        DATASTRUCT.('Spurious_FRET').File(i).Data.EA = DATASTRUCT.('Spurious_FRET').File(i).Data.YFPfret./DATASTRUCT.('Spurious_FRET').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
        DATASTRUCT.('Spurious_FRET').File(i).Data.ED = DATASTRUCT.('Spurious_FRET').File(i).Data.YFPfret./(DATASTRUCT.('Spurious_FRET').File(i).Data.YFPfret+ DATASTRUCT.('Spurious_FRET').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);


        ND = DATASTRUCT.('Spurious_FRET').File(i).Data.CFPdirect./(1-DATASTRUCT.('Spurious_FRET').File(i).Data.ED);
        NA = DATASTRUCT.('Spurious_FRET').File(i).Data.YFPdirect./(DATASTRUCT.PARAMETERS.gRatio*DATASTRUCT.PARAMETERS.fRatio);     
%         sum((ND./NA - DATASTRUCT.('Spurious_FRET').File(i).Data.EA./DATASTRUCT.('Spurious_FRET').File(i).Data.ED).^2)
       

        axes(MyHandles.Ax(4));
        MASK2 = MASK & (ND./NA<DATASTRUCT.MASKPARAMS.NDNA) & (NA./ND<DATASTRUCT.MASKPARAMS.NAND);
        line(ND(MASK2),DATASTRUCT.('Spurious_FRET').File(i).Data.EA(MASK2),'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor',MyColors(i,:));
        EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
        line([0 3E5],EASlope*[0 1E5],'linestyle','-','marker','none','linewidth',1,'color',[ 0 0 0]);
        if strcmp(get(MyHandles.CM(4),'checked'),'off')
            xlim([0 3E4]); ylim([-0.1 0.5])
        end
        axes(MyHandles.Ax(5));
        line(YFPTEMP(MASK2),DATASTRUCT.('Spurious_FRET').File(i).Data.ED(MASK2),'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor',MyColors(i,:));
        EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
        line([0 1E5],EDSlope*[0 1E5],'linestyle','-','marker','none','linewidth',1,'color',[ 0 0 0]);
        if strcmp(get(MyHandles.CM(5),'checked'),'off')
            xlim([0 3E4]); ylim([-0.1 0.5])
        end
    end    
    
% 
%     set(MyHandles.ParameterPanel.EASlope.Edit,'string',num2str(DATASTRUCT.PARAMETERS.EASlope));
%     set(MyHandles.ParameterPanel.EDSlope.Edit,'string',num2str(DATASTRUCT.PARAMETERS.EDSlope));
%     set(gcf,'userdata',DATASTRUCT);
    
function SetupConstructs()
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    SUBDATA = DATASTRUCT.('Constructs').File;
    
    sel = get(MyHandles.Constructs.FileList,'value');

    MyColors = jet(length(SUBDATA));
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;
    
    Gx = str2num(get(MyHandles.ParameterPanel.G.Edit,'string'));
    Fx = str2num(get(MyHandles.ParameterPanel.F.Edit,'string'));
    DATASTRUCT.PARAMETERS.fRatio = Fx;
    DATASTRUCT.PARAMETERS.gRatio = Gx;

    EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
    EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
    if isempty(EASlope) 
        EASlope = 0;
    end
    if isempty(EDSlope) 
        EDSlope = 0;
    end
        
    for i = sel
        CFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.CFP] - PARAMS.BG.CFP;
        YFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.YFP] - PARAMS.BG.YFP;
        FRETTEMP = [DATASTRUCT.('Constructs').File(i).Data.FRET] - PARAMS.BG.FRET;
        
        MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
        MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
        
        Fin = [CFPTEMP YFPTEMP FRETTEMP]';
        Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';
        
        DATASTRUCT.('Constructs').File(i).Data.CFPdirect  = Sout(:,1);
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect  = Sout(:,2);
        DATASTRUCT.('Constructs').File(i).Data.YFPfret = Sout(:,3);

        DATASTRUCT.('Constructs').File(i).Data.CFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPfret(MASK<1) = 0;
    end

    for i = sel
                  
        Kd = DATASTRUCT.('Constructs').File(i).FitParams.Kdeff;% str2num(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));
        nDnA = DATASTRUCT.('Constructs').File(i).FitParams.nDnA;%str2num(get(MyHandles.ParameterPanel.nDnA.Edit,'string'));    
        EAmax = DATASTRUCT.('Constructs').File(i).FitParams.EAmax;%str2num(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));
        EDmax = DATASTRUCT.('Constructs').File(i).FitParams.EDmax; %str2num(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));
        
        if length(sel)==1
            set(MyHandles.ParameterPanel.Kdeff.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.Kdeff)
            set(MyHandles.ParameterPanel.nDnA.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.nDnA)
            set(MyHandles.ParameterPanel.EAmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EAmax) 
            set(MyHandles.ParameterPanel.EDmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EDmax)
        end
        
        DATASTRUCT.('Constructs').File(i).Data.EA = DATASTRUCT.('Constructs').File(i).Data.YFPfret./DATASTRUCT.('Constructs').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
        DATASTRUCT.('Constructs').File(i).Data.ED = DATASTRUCT.('Constructs').File(i).Data.YFPfret./(DATASTRUCT.('Constructs').File(i).Data.YFPfret+ DATASTRUCT.('Constructs').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);


        
        ND = DATASTRUCT.('Constructs').File(i).Data.CFPdirect./(1-DATASTRUCT.('Constructs').File(i).Data.ED);
        NA = DATASTRUCT.('Constructs').File(i).Data.YFPdirect./(DATASTRUCT.PARAMETERS.gRatio*DATASTRUCT.PARAMETERS.fRatio);        
        MASK2 = MASK & (ND./NA<DATASTRUCT.MASKPARAMS.NDNA) & (NA./ND<DATASTRUCT.MASKPARAMS.NAND);
        
        EA_corrected = DATASTRUCT.('Constructs').File(i).Data.EA - ND*EASlope;
        ED_corrected = DATASTRUCT.('Constructs').File(i).Data.ED - NA*EDSlope; 
        
        % Compute Dfree / Afree
        
        Db = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*ND);
        Ab = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*NA);
        
%         Dfree = ((ND - Kd - NA)+sqrt((ND-Kd-NA).^2+4.*Kd.*NA))/2;
%         Afree = ((NA - Kd - ND)+sqrt((NA-Kd-ND).^2+4.*Kd.*ND))/2;    
        Dfree = ND.*(1-Db);
        Afree = NA.*(1-Ab);
        
        Dfree_F = Dfree(MASK2);
        Afree_F = Afree(MASK2);
        EA_corrected_F = EA_corrected(MASK2);
        ED_corrected_F = ED_corrected(MASK2);
        
        
        % binning data: 
        dB = 0.03;
        Fbcenter = dB/2:dB:1+dB;

        DCenter = Kd*Fbcenter./(1-Fbcenter);
        [DfreeCenter EACenter XSEM YSEM] = binData(Dfree_F, EA_corrected_F, DCenter);
        [AfreeCenter EDCenter XSEM YSEM] = binData(Afree_F, ED_corrected_F, DCenter);

       
        axes(MyHandles.Ax(4));
        cla        
            line(Dfree_F,EA_corrected_F,'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0]);
            try
                line(DfreeCenter, EACenter,'Marker','o','LineStyle','none','MarkerSize',5,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0.5 0.5]);
                
            end
            if strcmp(get(MyHandles.CM(4),'checked'),'off')
                ylim([-0.1 0.8]);xlim([0 1.2E5])
            end
            if ~isempty(EAmax)
                DfreeFit = logspace(-5, 7, 10000);
                EAFit = EAmax.*(DfreeFit./(DfreeFit+Kd));
                line(DfreeFit, EAFit, 'marker','none','linestyle','-','color',[1 0 0],'linewidth',2);
                line(DfreeFit, 0*DfreeFit,'marker','none','linestyle','-','color',[0.5 0.5 0.5],'linewidth',0.5);
                if strcmp(get(MyHandles.CM(4),'checked'),'off')
                    if EAmax > 0
                        ylim([-0.05 2*EAmax]);xlim([0 1.2E5])
                    else
                        ylim([-0.05 0.7]);xlim([0 1.2E5])     
                    end
                end
            end
%              set(gca,'Xscale','log')
        axes(MyHandles.Ax(5));
        cla
           line(Afree_F,ED_corrected_F,'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0]);
           if strcmp(get(MyHandles.CM(5),'checked'),'off')
                ylim([-0.05 3E5]);xlim([0 12E4])
           end
           try
            line(AfreeCenter, EDCenter,'Marker','o','LineStyle','none','MarkerSize',5,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 .5 .5]);
            end
            if ~isempty(EAmax)
                AfreeFit = logspace(-3, 7, 1000);
                EDFit = EDmax.*(AfreeFit./(AfreeFit+Kd/nDnA));
                line(AfreeFit, EDFit, 'marker','none','linestyle','-','color',[1 0 0],'linewidth',2);
                line(AfreeFit, 0*AfreeFit,'marker','none','linestyle','-','color',[0.5 0.5 0.5],'linewidth',0.5);
                if strcmp(get(MyHandles.CM(5),'checked'),'off')
                    if EDmax > 0
                        ylim([-0.05 2*EDmax]);xlim([0 1.2E5])
                    else
                        ylim([-0.05 0.7]);xlim([0 1.2E5])    
                    end
                end
            end
%              set(gca,'Xscale','log')
            
%             clipboard('copy',sprintf('%7.6f\t%7.6f',sum(TEMPGATE>rx)/length(TEMPGATE),trimmean(NANDRATIO(MASK_X>0),95)))
             

    end        
    
    
function Fout = TransformFluorescence(Fin,RA1,RD1,RD2)
    T = [RD1 0 0; -RA1*RD2 RA1 0; (RA1*RD2-RD1) -RA1 1];
    Fout = T*Fin;
    
function varargout = removeNaN(varargin)
    for i = 1:nargin
        MASK(:,i) = isnan(varargin{i});        
    end
    MASK = sum(MASK,2)<1;
    for i = 1:nargin
        varargout{i} = varargin{i}(MASK);
    end
    
function UpdateCalib(tagIndex)
%     DATASTRUCT = get(gcf,'userdata');  
    global DATASTRUCT;
    MyHandles = DATASTRUCT.Handles;
    switch tagIndex
        case 1
            DATASTRUCT.PARAMETERS.BG.CFP = str2num(get(MyHandles.ParameterPanel.BGValsC.Edit,'string'));
            if ~DATASTRUCT.WORKFLOW.bMaskSet
                DATASTRUCT.MASKPARAMS.CFP.Min = DATASTRUCT.PARAMETERS.BG.CFP_CUTOFF;
            end
            UpdatePlots('Background');    
        case 2
            DATASTRUCT.PARAMETERS.BG.YFP = str2num(get(MyHandles.ParameterPanel.BGValsY.Edit,'string')); 
            if ~DATASTRUCT.WORKFLOW.bMaskSet
                DATASTRUCT.MASKPARAMS.YFP.Min = DATASTRUCT.PARAMETERS.BG.YFP_CUTOFF;
            end
            UpdatePlots('Background');
        case 3
            DATASTRUCT.PARAMETERS.BG.FRET = str2num(get(MyHandles.ParameterPanel.BGValsF.Edit,'string'));
            if ~DATASTRUCT.WORKFLOW.bMaskSet
                DATASTRUCT.MASKPARAMS.FRET.Min = DATASTRUCT.PARAMETERS.BG.FRET_CUTOFF;
            end
            UpdatePlots('Background');
        case 4
            DATASTRUCT.PARAMETERS.RA1 = str2num(get(MyHandles.ParameterPanel.RA1.Edit,'string'));
            UpdatePlots('Acceptor');
        case 5
            DATASTRUCT.PARAMETERS.RD1 = str2num(get(MyHandles.ParameterPanel.RD1.Edit,'string'));       
            UpdatePlots('Donor');
        case 6
            DATASTRUCT.PARAMETERS.RD2 = str2num(get(MyHandles.ParameterPanel.RD2.Edit,'string'));    
            UpdatePlots('Donor');
        case 8
            DATASTRUCT.PARAMETERS.fRatio = str2num(get(MyHandles.ParameterPanel.F.Edit,'string')); 
 
        case 7
            DATASTRUCT.PARAMETERS.gRatio = str2num(get(MyHandles.ParameterPanel.G.Edit,'string'));      
        % copy this 
        case 9
            DATASTRUCT.PARAMETERS.EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));  
            DATASTRUCT.PARAMETERS.EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));  
            UpdatePlots('Spurious_FRET');
        case 10
            DATASTRUCT.PARAMETERS.EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));  
            DATASTRUCT.PARAMETERS.EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));  
            UpdatePlots('Spurious_FRET');
            
         % copy this     
        case 11
            % 
            sel = get(MyHandles.Constructs.FileList,'value');
            DATASTRUCT.Constructs.File(sel).FitParams.Kdeff = str2num(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));
%             set(gcf,'userdata',DATASTRUCT);  
            UpdatePlots('Constructs');
        case 14
            sel = get(MyHandles.Constructs.FileList,'value');
            DATASTRUCT.Constructs.File(sel).FitParams.nDnA = str2num(get(MyHandles.ParameterPanel.nDnA.Edit,'string'));
%             set(gcf,'userdata',DATASTRUCT);  
            UpdatePlots('Constructs');
            
        case 12
            sel = get(MyHandles.Constructs.FileList,'value');
            DATASTRUCT.Constructs.File(sel).FitParams.EAmax = str2num(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));
%             set(gcf,'userdata',DATASTRUCT);  
            UpdatePlots('Constructs');
            
        case 13
            sel = get(MyHandles.Constructs.FileList,'value');
            DATASTRUCT.Constructs.File(sel).FitParams.EDmax = str2num(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));
%             set(gcf,'userdata',DATASTRUCT);  
            UpdatePlots('Constructs');

            
            
    end
%     set(gcf,'userdata',DATASTRUCT);  
        
function PARAMS = getPARAMETERS(MyHandles)
    PARAMS.BG.CFP   = str2double(get(MyHandles.ParameterPanel.BGValsC.Edit,'string'));
    PARAMS.BG.YFP   = str2double(get(MyHandles.ParameterPanel.BGValsY.Edit,'string'));  
    PARAMS.BG.FRET  = str2double(get(MyHandles.ParameterPanel.BGValsF.Edit,'string'));
    
    PARAMS.BG.CFP_CUTOFF   = str2double(get(MyHandles.ParameterPanel.BGValsC_cut.Edit,'string'));
    PARAMS.BG.YFP_CUTOFF   = str2double(get(MyHandles.ParameterPanel.BGValsY_cut.Edit,'string'));  
    PARAMS.BG.FRET_CUTOFF  = str2double(get(MyHandles.ParameterPanel.BGValsF_cut.Edit,'string'));
    
    
    PARAMS.RA1      = str2double(get(MyHandles.ParameterPanel.RA1.Edit,'string'));
    PARAMS.RD1      = str2double(get(MyHandles.ParameterPanel.RD1.Edit,'string'));
    PARAMS.RD2      = str2double(get(MyHandles.ParameterPanel.RD2.Edit,'string'));
    PARAMS.fRatio   = str2double(get(MyHandles.ParameterPanel.F.Edit,'string'));
    PARAMS.gRatio   = str2double(get(MyHandles.ParameterPanel.G.Edit,'string'));  
    
    PARAMS.EASlope  = str2double(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));  
    PARAMS.EDSlope  = str2double(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));  
    
    PARAMS.Kdeff    = str2double(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));  
    PARAMS.EAmax    = str2double(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));  
    PARAMS.EDmax    = str2double(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));  
    
    


function setPARAMETERS(PARAMS, MyHandles)
%     PARAMS = getPARAMETERS(MyHandles);
    set(MyHandles.ParameterPanel.BGValsC.Edit,'string',num2str(PARAMS.BG.CFP));
    set(MyHandles.ParameterPanel.BGValsY.Edit,'string',num2str(PARAMS.BG.YFP));
    set(MyHandles.ParameterPanel.BGValsF.Edit,'string',num2str(PARAMS.BG.FRET));
    
    set(MyHandles.ParameterPanel.BGValsC_cut.Edit,'string',num2str(PARAMS.BG.CFP_CUTOFF));
    set(MyHandles.ParameterPanel.BGValsY_cut.Edit,'string',num2str(PARAMS.BG.YFP_CUTOFF));
    set(MyHandles.ParameterPanel.BGValsF_cut.Edit,'string',num2str(PARAMS.BG.FRET_CUTOFF));
    
    set(MyHandles.ParameterPanel.RA1.Edit,'string',num2str(PARAMS.RA1));
    set(MyHandles.ParameterPanel.RD1.Edit,'string',num2str(PARAMS.RD1));
    set(MyHandles.ParameterPanel.RD2.Edit,'string',num2str(PARAMS.RD2));
        
    set(MyHandles.ParameterPanel.F.Edit,'string',num2str(PARAMS.fRatio));
    set(MyHandles.ParameterPanel.G.Edit,'string',num2str(PARAMS.gRatio)); 
    
    set(MyHandles.ParameterPanel.EASlope.Edit,'string',num2str(PARAMS.EASlope));
    set(MyHandles.ParameterPanel.EDSlope.Edit,'string',num2str(PARAMS.EDSlope)); 
    
    
function [slope] = iterativeScalarFit(XX,YY)
    try
        g = fittype('a*u+b','independent','u');
        [~,~,fitinfo]=fit(XX,YY,g,'StartPoint',[1 1]);
        residuals = fitinfo.residuals;
        I = abs( residuals) > 1.5 * std( residuals );
        outliers = excludedata(XX,YY,'indices',I);
        fit2 = fit(XX,YY,g,'StartPoint',[1 1],'Exclude',outliers);
        slope = fit2.a;
        intercept = fit2.b;
    catch
       p = polyfit(XX,YY,1);
       RA1=p(1); 
    end
        
function [Gx Fx] = getFRETStd(SUBDATA, MainFig, FigNew)
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
        
    Ax = axes('units','normalized','position', [0.2 0.2 0.6 0.6]);
    
    nFiles = length(SUBDATA.File);
    
    Gx = str2num(get(MyHandles.ParameterPanel.G.Edit,'string'));
    Fx = str2num(get(MyHandles.ParameterPanel.F.Edit,'string'));
    MyColors = jet(nFiles);
    
    SubG.Label = uicontrol('parent',FigNew, 'style','text', 'units','normalized','position',[0.1 0.1 .1 .05],'string','slope','BackgroundColor',[1 1 1]);
    SubG.Edit = uicontrol('parent',FigNew, 'style','edit', 'units','normalized','position',[0.2 0.1 .1 .05],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateFGPlot'',' num2str(FigNew) ')']);
    SubF.Label = uicontrol('parent',FigNew, 'style','text', 'units','normalized','position',[0.4 0.1 .1 .05],'string','intercept','BackgroundColor',[1 1 1]);
    SubF.Edit = uicontrol('parent',FigNew, 'style','edit', 'units','normalized','position',[0.5 0.1 .1 .05],'string','','BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateFGPlot'',' num2str(FigNew) ')']);
    bReplace = 1;
    if ~isempty(Gx) && ~isempty(Fx) && ~(Gx==0) && ~(Fx==0)  
        button = questdlg('Replace Calibration Data?', 'FACS FRET', 'Yes', 'No', 'No');
        switch button
            case 'Yes'
                bReplace = 1;
            otherwise
                bReplace = 0;
        end
    else
        bReplace = 1;
    end     
    CALIBDATA_median = zeros(length(SUBDATA.File),2);
    for i = 1:nFiles
        Rx = SUBDATA.File(i).Data.YFPdirect./SUBDATA.File(i).Data.CFPdirect;
        Ry = SUBDATA.File(i).Data.YFPfret./SUBDATA.File(i).Data.CFPdirect;
        MASKx = SUBDATA.File(i).Data.CFPdirect==0;
        Rx = Rx(~MASKx);
        Ry = Ry(~MASKx);
        [Rx Ry] = removeNaN(Rx,Ry);
        CALIBDATA_median(i,:) = [median(Rx), median(Ry)];
        line(Rx, Ry, 'Marker','o','LineStyle','none','MarkerSize',1,'MarkerEdgeColor','none','MarkerFaceColor',MyColors(i,:));
    end
    line(CALIBDATA_median(:,1),CALIBDATA_median(:,2),'Marker','o','LineStyle','none','MarkerSize',5,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0]);    
    xlim([0 0.2]);ylim([-2 5]);



    if bReplace>0  
        b = robustfit(CALIBDATA_median(:,1),CALIBDATA_median(:,2),'fair',1.4);
        Rx_Fit = linspace(0,0.2,2);
        Ry_Fit = b(1) + b(2)*Rx_Fit;           
        Gx = 1/b(2);
        Fx = -b(1);
        set(SubG.Edit, 'string', b(2));
        set(SubF.Edit, 'string', b(1));     
    else
        Rx_Fit = linspace(0,0.2,2);                     
        slope = 1/Gx;
        intercept = -Fx;        
        Ry_Fit = intercept + slope*Rx_Fit;
        set(SubG.Edit, 'string', slope);
        set(SubF.Edit, 'string', intercept);    
        
    end
    AFit = line(Rx_Fit, Ry_Fit, 'linestyle','-','linewidth',2,'Marker','none');
    set(FigNew, 'userdata',{Ax AFit MainFig SubG.Edit SubF.Edit});
    UpdateFGPlot(FigNew)
   
    
function UpdateFGPlot(FigNew)
    RelevantHandles = get(FigNew,'userdata');
    Ax = RelevantHandles{1};
    HandleLine = RelevantHandles{2};
    MainFig = RelevantHandles{3};
    HandleG = RelevantHandles{4};
    HandleF = RelevantHandles{5};
    
    slope = str2num(get(HandleG,'string'));
    Gx = 1/slope;
    intercept  = str2num(get(HandleF,'string'));    
    Fx = -intercept;
    Rx = get(HandleLine,'XData');
    Ry = -Fx + 1/Gx*Rx;
    set(HandleLine,'YData',Ry);
    
    global DATASTRUCT
%     DATASTRUCT = get(MainFig,'userData');    
    MyHandles = DATASTRUCT.Handles;
    set(MyHandles.ParameterPanel.G.Edit,'string',num2str(Gx));
    set(MyHandles.ParameterPanel.F.Edit,'string',num2str(Fx));
    
    
function [XCenter YCenter XSEM YSEM] = binData(X, Y, XBin)

NAN_MASK = isnan(X) | isnan(Y) | (X<0) |(Y<-0.04);
X = X(~NAN_MASK);
Y = Y(~NAN_MASK);

XCenter = zeros(length(XBin),1);
YCenter = zeros(length(XBin),1);
YSEM = XCenter;
XSEM = XCenter;
cutoff = 99;
MASK = X<XBin(1);
XCenter(1) = median(X(MASK));
YCenter(1) = median(Y(MASK));
XSEM(1) = std(X(MASK))./sqrt(length(X(MASK)));
YSEM(1) = std(Y(MASK))./sqrt(length(Y(MASK)));

for i = 2:length(XBin)
    MASK = (X<XBin(i))&(X>=XBin(i-1));
    XCenter(i) = median(X(MASK));
    YCenter(i) = median(Y(MASK));
    XSEM(i) = std(X(MASK))./sqrt(length(X(MASK)));
    YSEM(i) = std(Y(MASK))./sqrt(length(Y(MASK)));
end

function Filtering()
global DATASTRUCT
MyHandles = DATASTRUCT.Handles;
PARAMS = getPARAMETERS(MyHandles);
if isfield(DATASTRUCT.Handles,'MASKFig') && ~isempty(DATASTRUCT.Handles.MASKFig)
   figure(DATASTRUCT.Handles.MASKFig)
else
   DATASTRUCT.Handles.MASKFig = figure('name',['MASK Parameters'], 'NumberTitle','off','units','normalized', 'position', [0 0 0.2 0.4],'color', [1 1 1],'toolbar','none','menubar','none','DeleteFcn',['FACS_FRET(''CloseMASKWindow'')']);

end

if ~DATASTRUCT.WORKFLOW.bMaskSet
    DATASTRUCT.MASKPARAMS.CFP.Min = PARAMS.BG.CFP_CUTOFF;
    DATASTRUCT.MASKPARAMS.CFP.Max = 2.5E5;
    DATASTRUCT.MASKPARAMS.YFP.Min = PARAMS.BG.YFP_CUTOFF;
    DATASTRUCT.MASKPARAMS.YFP.Max = 2.5E5;
    DATASTRUCT.MASKPARAMS.FRET.Min = PARAMS.BG.FRET_CUTOFF;
    DATASTRUCT.MASKPARAMS.FRET.Max = 2.5E5;
    DATASTRUCT.MASKPARAMS.NDNA = 50;
    DATASTRUCT.MASKPARAMS.NAND = 50;
end

BGW = 1/5;BGS = BGW/15; j =2; BGH = 1/12;
FilterParams.CFP.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*1         1-BGH-(BGH+BGS)*(j-1) BGW/2 BGH],'string','CFP :','BackgroundColor',[0.8 1 1],'HorizontalAlignment','left');
FilterParams.CFPMin.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*2+BGW*0.5 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string','Minimum','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
FilterParams.CFPMin.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*3+BGW*1.5 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string',num2str(DATASTRUCT.MASKPARAMS.CFP.Min),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',1)']);
FilterParams.CFPMax.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*4+2.5*BGW 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string','Maximum','BackgroundColor',[1 1 1],'HorizontalAlignment','center'); 
FilterParams.CFPMax.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*5+3.5*BGW 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string',num2str(DATASTRUCT.MASKPARAMS.CFP.Max),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',2)']);

j = j+1;
FilterParams.YFP.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*1         1-BGH-(BGH+BGS)*(j-1) BGW/2 BGH],'string','YFP :','BackgroundColor',[1 1 0.8],'HorizontalAlignment','left');
FilterParams.YFPMin.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*2+BGW*0.5 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string','Minimum','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
FilterParams.YFPMin.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*3+BGW*1.5 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string',num2str(DATASTRUCT.MASKPARAMS.YFP.Min),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',3)']);
FilterParams.YFPMax.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*4+2.5*BGW 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string','Maximum','BackgroundColor',[1 1 1],'HorizontalAlignment','center'); 
FilterParams.YFPMax.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*5+3.5*BGW 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string',num2str(DATASTRUCT.MASKPARAMS.YFP.Max),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',4)']);

j = j+1;
FilterParams.FRET.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*1         1-BGH-(BGH+BGS)*(j-1) BGW/2 BGH],'string','FRET :','BackgroundColor',[0.8 1 0.8],'HorizontalAlignment','left');
FilterParams.FRETMin.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*2+BGW*0.5 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string','Minimum','BackgroundColor',[1 1 1],'HorizontalAlignment','center');
FilterParams.FRETMin.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*3+BGW*1.5 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string',num2str(DATASTRUCT.MASKPARAMS.FRET.Min),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',5)']);
FilterParams.FRETMax.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS*4+2.5*BGW 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string','Maximum','BackgroundColor',[1 1 1],'HorizontalAlignment','center'); 
FilterParams.FRETMax.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*5+3.5*BGW 1-BGH-(BGH+BGS)*(j-1) BGW   BGH],'string',num2str(DATASTRUCT.MASKPARAMS.FRET.Max),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',6)']);

j = j+1;

FilterParams.NAND.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS 1-BGH-(BGH+BGS)*(j-1) BGW BGH],'string','NA/ND Maximum','BackgroundColor',[1 1 1],'HorizontalAlignment','left'); 
FilterParams.NAND.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*1+BGW 1-BGH-(BGH+BGS)*(j-1) BGW BGH],'string',num2str(DATASTRUCT.MASKPARAMS.NAND),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',7)']);

j = j+1;
FilterParams.NDNA.Label = uicontrol( 'style','text', 'units','normalized','position',[BGS 1-BGH-(BGH+BGS)*(j-1) BGW BGH],'string','ND/NA Maximum','BackgroundColor',[1 1 1],'HorizontalAlignment','left'); 
FilterParams.NDNA.Edit = uicontrol('style','edit', 'units','normalized','position',  [BGS*1+BGW 1-BGH-(BGH+BGS)*(j-1) BGW BGH],'string',num2str(DATASTRUCT.MASKPARAMS.NDNA),'BackgroundColor',[1 1 1],'CallBack',['FACS_FRET(''UpdateMask'',8)']);

j = j+1;
DATASTRUCT.Handles.FilterParams = FilterParams;
% FilterParams.ApplyMask = uicontrol('style','pushbutton', 'units','normalized','position',[BGS 1-BGH-(BGH+BGS)*(j-1) BGW BGH],'string','Apply Mask','BackgroundColor',[.8 1 .8], 'CallBack',['FACS_FRET(''ApplyMask'')'],'SelectionHighlight','off');

function UpdateMask(tagIndex)
    global DATASTRUCT;
    MyHandles = DATASTRUCT.Handles;
    switch tagIndex
        case 1
            DATASTRUCT.MASKPARAMS.CFP.Min = str2num(get(DATASTRUCT.Handles.FilterParams.CFPMin.Edit,'string'));
    
        case 2
            DATASTRUCT.MASKPARAMS.CFP.Max = str2num(get(DATASTRUCT.Handles.FilterParams.CFPMax.Edit,'string'));
   
        case 3
             DATASTRUCT.MASKPARAMS.YFP.Min = str2num(get(DATASTRUCT.Handles.FilterParams.YFPMin.Edit,'string'));
    
        case 4
            DATASTRUCT.MASKPARAMS.YFP.Max = str2num(get(DATASTRUCT.Handles.FilterParams.YFPMax.Edit,'string'));
   
        case 5
             DATASTRUCT.MASKPARAMS.FRET.Min = str2num(get(DATASTRUCT.Handles.FilterParams.FRETMin.Edit,'string'));
        case 6     
            DATASTRUCT.MASKPARAMS.FRET.Max = str2num(get(DATASTRUCT.Handles.FilterParams.FRETMax.Edit,'string'));
        case 7
            DATASTRUCT.MASKPARAMS.NAND = str2num(get(DATASTRUCT.Handles.FilterParams.NAND.Edit,'string'));
        case 8
            DATASTRUCT.MASKPARAMS.NDNA = str2num(get(DATASTRUCT.Handles.FilterParams.NDNA.Edit,'string'));
    end
    DATASTRUCT.WORKFLOW.bMaskSet = 1;
    UpdatePlots('')
    figure(DATASTRUCT.Handles.MASKFig)
    
function CloseMASKWindow()
global DATASTRUCT
DATASTRUCT.Handles.MASKFig = [];
        

function ShowHideRawData()
% MainFig = gcf;
% DATASTRUCT = get(MainFig,'userdata');
global DATASTRUCT
UpdatePlots('')
PanelAx2Pos1 = [0.2200    0.0100    0.7700    0.4550];
PanelAx2Pos2 = [0.2200    0.0100    0.7700    0.73];
if get(DATASTRUCT.Handles.AnalysisPanel.ShowHidePanel,'value')>0
    set(DATASTRUCT.Handles.PanelAx,'visible','on')
    set(DATASTRUCT.Handles.PanelAx2,'Position',PanelAx2Pos1);
    set(DATASTRUCT.Handles.AnalysisPanel.ShowHidePanel,'BackgroundColor',[0.7 1 0.7])
else
    set(DATASTRUCT.Handles.PanelAx,'visible','off')
    set(DATASTRUCT.Handles.PanelAx2,'Position',PanelAx2Pos2);
    set(DATASTRUCT.Handles.AnalysisPanel.ShowHidePanel,'BackgroundColor',[1 0.7 0.7])
end

function ClearAllAxes()
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;

    for i = 1:length(MyHandles.Ax)
        cla(MyHandles.Ax(i));
    end


function varargout = filterNaN(varargin)
MaskNaN = zeros(size(varargin{1},1),length(varargin));
try
    for i = 1: length(varargin)
        TempMask = MaskNaN | isnan(varargin{i});
    end
    for i  = 1:length(vararin)
        Temp = varargin{i};
        varargout{i} = Temp(~MaskNaN);
    end
catch EM
    disp(EM);
    varargout = {};
    return
end

function ChangeAxLimits(varargin)
    global DATASTRUCT
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    MyHandles = DATASTRUCT.Handles;
    
    if strcmp(get(MyHandles.CM(varargin{1}),'checked'),'off')
        XLIMIT = get(MyHandles.Ax(varargin{1}),'xlim');
        YLIMIT = get(MyHandles.Ax(varargin{1}),'ylim');    
        XLABEL = get(get(MyHandles.Ax(varargin{1}),'xlabel'),'string');
        if isempty(XLABEL)
            XLABEL = 'X-axis';
        end
        YLABEL = get(get(MyHandles.Ax(varargin{1}),'ylabel'),'string');
        if isempty(XLABEL)
            XLABEL = 'Y-axis';
        end
        prompt={['Define ' XLABEL ' limits:'],...
            ['Define ' YLABEL ' limits:']};

        name='Setup Axes';
        numlines=1;
        defaultanswer={num2str(round(XLIMIT*10)/10),num2str(round(YLIMIT*10)/10)};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        if ~isempty(answer)
            NEWXLIMIT = str2num(answer{1});
            NEWYLIMIT = str2num(answer{2});
            set(MyHandles.CM(varargin{1}),'checked','on')
        else
            NEWXLIMIT = XLIMIT;
            NEWYLIMIT = YLIMIT;   
        end
        set(MyHandles.Ax(varargin{1}),'xlim',NEWXLIMIT,'ylim',NEWYLIMIT)
    else
        set(MyHandles.CM(varargin{1}),'checked','off')
        UpdatePlots('');
    end
    
    
    
function SaveCalibration(varargin)
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    PARAMS = getPARAMETERS(MyHandles);
    [fname path] = uiputfile('*.par', 'Save Parameters',DATASTRUCT.ActivePath);
    
    
    
    save([path fname],'-struct','PARAMS');
    
        
function LoadCalibration(varargin)
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    [fname path] = uigetfile('*.par', 'Save Parameters',DATASTRUCT.ActivePath);
    if ~fname==0
        PARAMS = load([path fname],'-mat');
        setPARAMETERS(PARAMS, MyHandles)
        UpdatePlots('')
    end
    
function ExportFigure(varargin)
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    
    SUBDATA = DATASTRUCT.('Constructs').File;
    
    sel = get(MyHandles.Constructs.FileList,'value');
    
    MyColors = jet(length(SUBDATA));
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;
    Gx = PARAMS.gRatio;
    Fx = PARAMS.fRatio;
    DATASTRUCT.PARAMETERS.fRatio = Fx;
    DATASTRUCT.PARAMETERS.gRatio = Gx;
    
    EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
    EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
    if isempty(EASlope) 
        EASlope = 0;
    end
    if isempty(EDSlope) 
        EDSlope = 0;
    end
        
    for i = sel
        CFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.CFP] - PARAMS.BG.CFP;
        YFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.YFP] - PARAMS.BG.YFP;
        FRETTEMP = [DATASTRUCT.('Constructs').File(i).Data.FRET] - PARAMS.BG.FRET;
        
        MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
        MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
            
        Fin = [CFPTEMP YFPTEMP FRETTEMP]';
        Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';
        
        DATASTRUCT.('Constructs').File(i).Data.CFPdirect  = Sout(:,1);
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect  = Sout(:,2);
        DATASTRUCT.('Constructs').File(i).Data.YFPfret = Sout(:,3);

        DATASTRUCT.('Constructs').File(i).Data.CFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPfret(MASK<1) = 0;
    end

    for i = sel
                  
        Kd = DATASTRUCT.('Constructs').File(i).FitParams.Kdeff;% str2num(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));
        nDnA = DATASTRUCT.('Constructs').File(i).FitParams.nDnA;%str2num(get(MyHandles.ParameterPanel.nDnA.Edit,'string'));    
        EAmax = DATASTRUCT.('Constructs').File(i).FitParams.EAmax;%str2num(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));
        EDmax = DATASTRUCT.('Constructs').File(i).FitParams.EDmax; %str2num(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));
        
        if length(sel)==1
            set(MyHandles.ParameterPanel.Kdeff.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.Kdeff)
            set(MyHandles.ParameterPanel.nDnA.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.nDnA)
            set(MyHandles.ParameterPanel.EAmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EAmax) 
            set(MyHandles.ParameterPanel.EDmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EDmax)
        end
        
        DATASTRUCT.('Constructs').File(i).Data.EA = DATASTRUCT.('Constructs').File(i).Data.YFPfret./DATASTRUCT.('Constructs').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
        DATASTRUCT.('Constructs').File(i).Data.ED = DATASTRUCT.('Constructs').File(i).Data.YFPfret./(DATASTRUCT.('Constructs').File(i).Data.YFPfret+ DATASTRUCT.('Constructs').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);


        
        ND = DATASTRUCT.('Constructs').File(i).Data.CFPdirect./(1-DATASTRUCT.('Constructs').File(i).Data.ED);
        NA = DATASTRUCT.('Constructs').File(i).Data.YFPdirect./(DATASTRUCT.PARAMETERS.gRatio*DATASTRUCT.PARAMETERS.fRatio);        
        MASK2 = MASK & (ND./NA<DATASTRUCT.MASKPARAMS.NDNA) & (NA./ND<DATASTRUCT.MASKPARAMS.NAND);
        
        EA_corrected = DATASTRUCT.('Constructs').File(i).Data.EA - ND*EASlope;
        ED_corrected = DATASTRUCT.('Constructs').File(i).Data.ED - NA*EDSlope; 
        
        EA_corrected_F = EA_corrected(MASK2);
        ED_corrected_F = ED_corrected(MASK2);
        ND_F = ND(MASK2);
        NA_F = NA(MASK2);
        
        
        
        Db = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*ND);
        Ab = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*NA);
           
        Dfree = ND.*(1-Db);
        Afree = NA.*(1-Ab);
        
        Dfree_F = Dfree(MASK2);
        Afree_F = Afree(MASK2);
        
        % binning data: 
        dB = 0.03;
        Fbcenter = dB/2:dB:1+dB;

        DCenter = Kd*Fbcenter./(1-Fbcenter);
        [DfreeCenter EACenter XSEM_A YSEM_A] = binData(Dfree_F, EA_corrected_F, DCenter);
        [AfreeCenter EDCenter XSEM_D YSEM_D] = binData(Afree_F, ED_corrected_F, DCenter);

        
        XLIMIT = get(MyHandles.Ax(5),'xlim');
        YLIMIT = get(MyHandles.Ax(5),'ylim');
        answer = inputdlg({'X-axis Limit:','Y-axis Limit:'},'Axes Limits',1,{num2str(XLIMIT),num2str(YLIMIT)});
        XLIMIT = str2num(answer{1});YLIMIT = str2num(answer{2});
        
        figure; 
            Ax1 = axes('units','centimeters', 'position',[2 2 5 3],'TickDir','out','TickLen',[0.03 0.03]);  
            try
                line(DfreeCenter, EACenter,'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0.5 0.5]);
                % error bars on each bin
            end
            line(Dfree_F,EA_corrected_F,'Marker','o','LineStyle','none','MarkerSize',1,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0]);

            DfreeFit = [0 logspace(-5, 7, 10000)];
            EAFit = EAmax.*(DfreeFit./(DfreeFit+Kd));
            line(DfreeFit, EAFit, 'marker','none','linestyle','-','color',[1 0 0],'linewidth',0.5);
            line([0 DfreeFit(end)], [0 DfreeFit(end)]*0,'marker','none','linestyle','-','color',[0.5 0.5 0.5],'linewidth',0.5);           


            set(Ax1,'xlim',XLIMIT,'ylim',YLIMIT);
            xlabel('Dfree (a.u.)');ylabel('E_A')
                
            Ax2 = axes('units','centimeters', 'position',[9 2 5 3],'TickDir','out','TickLen',[0.03 0.03]);  
            try
                line(AfreeCenter, EDCenter,'Marker','o','LineStyle','none','MarkerSize',2,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 .5 .5]);
            end
            line(Afree_F,ED_corrected_F,'Marker','o','LineStyle','none','MarkerSize',1,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0]);

            AfreeFit = [0 logspace(-5, 7, 1000)];
            EDFit = EDmax.*(AfreeFit./(AfreeFit+Kd/nDnA));
            line(AfreeFit, EDFit, 'marker','none','linestyle','-','color',[1 0 0],'linewidth',0.5);
            line([0 AfreeFit(end)], [0 0],'marker','none','linestyle','-','color',[0.5 0.5 0.5],'linewidth',0.5); 
            set(Ax2,'xlim',XLIMIT,'ylim',YLIMIT);
            xlabel('Afree (a.u.)');ylabel('E_D')
            
    end        

    
    function [initFitParam CI] = AutoKDfit(ED,EA,ND,NA, initFitParam)
        X = [ND NA];
        Kd = initFitParam(1);
        EAmax = initFitParam(2);
        EDmax = initFitParam(3);
        
        button = questdlg({'Please choose parameters to optimize for Autofit:' '' '(1) FreeFit -- optimizes Kd, EAmax, and EDmax' ...
            '(2) Kd only -- optimizes Kd given fixed EAmax and EDmax values' '(3) EMax -- optimize EAmax and EDmax given fixed Kd' '' ''},'Fit Params','Free Fit','Kd Only','Emax Only','Free Fit');
        
        switch button
            case 'Free Fit'
                options = statset('Robust','on');
                [beta,resid,J,sigma]  = nlinfit(X,ED,@EDsinglebinding, [Kd,EDmax],options);
                cix = nlparci(beta,resid,'covar',sigma) ;
                [EAmax,resid2,J2,sigma2] = nlinfit(X,EA,@(phi,X)EAsinglebindingED(phi,beta(1),X),beta(2));
                cix2 = nlparci(EAmax,resid2,'covar',sigma2);
                initFitParam(1) = beta(1);
                initFitParam(2) = EAmax;
                initFitParam(3) = beta(2);
                CI(1,:) = cix(1,:);CI(2,:) = cix2(1,:); CI(3,:) = cix(2,:);
                
            case 'Kd Only'
                options = statset('Robust','off');
                [beta,resid,J,sigma]  = nlinfit(X,ED,@(phi,X)EDsinglebinding([phi,EDmax],X), Kd, options);
                cix = nlparci(beta,resid,'covar',sigma);
                initFitParam(1) = beta(1);
                initFitParam(2) = EAmax;
                initFitParam(3) = EDmax;
                CI(1,:) = cix(1,:);CI(2,:) = [EAmax EAmax]; CI(3,:) = [EDmax EDmax];
                
            case 'Emax Only'
                options = statset('Robust','off');
                [beta,resid,J,sigma]  = nlinfit(X,ED,@(phi,X)EDsinglebinding([Kd,phi],X), EDmax, options);
                cix = nlparci(beta,resid,'covar',sigma) ;
                [EAmax,resid2,J2,sigma2] = nlinfit(X,EA,@(phi,X)EAsinglebindingED(phi,Kd,X),beta);
                cix2 = nlparci(EAmax,resid2,'covar',sigma2);
                initFitParam(1) = Kd;
                initFitParam(2) = EAmax;
                initFitParam(3) = beta;
                CI(1,:) = [Kd Kd];CI(2,:) = cix2(1,:); CI(3,:) = cix(1,:);
            otherwise
                CI = [Kd Kd; EAmax EAmax; EDmax EDmax];
        end
        
    function EDfit = EDsinglebinding(beta,X)
        ND = X(:,1);NA = X(:,2);
        Kd = beta(1);
        EDmax = beta(2);
        
        % Compute Dfree / Afree
        Db = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*ND);
        Ab = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*NA);
        
        Dfree = ND.*(1-Db);
        Afree = NA.*(1-Ab);
% 
%         Dfree = ((ND - Kd - NA)+sqrt((ND-Kd-NA).^2+4.*Kd.*NA))/2;
%         Afree = ((NA - Kd - ND)+sqrt((NA-Kd-ND).^2+4.*Kd.*ND))/2;  
        EDfit = EDmax.*Afree./(Afree+Kd);
 
    function EAfit = EAsinglebindingED(EAmax,Kd,X)
        ND = X(:,1);NA = X(:,2);
        % Compute Dfree / Afree
        Db = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*ND);
        Ab = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*NA);
        
        Dfree = ND.*(1-Db);
        Afree = NA.*(1-Ab);
%         
%         Dfree = ((ND - Kd - NA)+sqrt((ND-Kd-NA).^2+4.*Kd.*NA))/2;
%         Afree = ((NA - Kd - ND)+sqrt((NA-Kd-ND).^2+4.*Kd.*ND))/2;    
        EAfit = EAmax.*Dfree./(Dfree+Kd);
        
    function AutoFitCaller(varargin)
%     MainFig = gcf;
    global DATASTRUCT
%     DATASTRUCT = get(MainFig,'userdata');
    MyHandles = DATASTRUCT.Handles;
    
    SUBDATA = DATASTRUCT.('Constructs').File;
    
    sel = get(MyHandles.Constructs.FileList,'value');
    
    MyColors = jet(length(SUBDATA));
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;
    Gx = PARAMS.gRatio;
    Fx = PARAMS.fRatio;
    DATASTRUCT.PARAMETERS.fRatio = Fx;
    DATASTRUCT.PARAMETERS.gRatio = Gx;
    
    EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
    EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
    if isempty(EASlope) 
        EASlope = 0;
    end
    if isempty(EDSlope) 
        EDSlope = 0;
    end
        
    for i = sel
        CFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.CFP] - PARAMS.BG.CFP;
        YFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.YFP] - PARAMS.BG.YFP;
        FRETTEMP = [DATASTRUCT.('Constructs').File(i).Data.FRET] - PARAMS.BG.FRET;
        
        Cutoff = 1;
        MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
        MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
            
        Fin = [CFPTEMP YFPTEMP FRETTEMP]';
        Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';
        
        DATASTRUCT.('Constructs').File(i).Data.CFPdirect  = Sout(:,1);
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect  = Sout(:,2);
        DATASTRUCT.('Constructs').File(i).Data.YFPfret = Sout(:,3);

        DATASTRUCT.('Constructs').File(i).Data.CFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPfret(MASK<1) = 0;
    end

    for i = sel
        Kd = DATASTRUCT.('Constructs').File(i).FitParams.Kdeff;% str2num(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));
        nDnA = DATASTRUCT.('Constructs').File(i).FitParams.nDnA;%str2num(get(MyHandles.ParameterPanel.nDnA.Edit,'string'));    
        EAmax = DATASTRUCT.('Constructs').File(i).FitParams.EAmax;%str2num(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));
        EDmax = DATASTRUCT.('Constructs').File(i).FitParams.EDmax; %str2num(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));
        
        if length(sel)==1
            set(MyHandles.ParameterPanel.Kdeff.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.Kdeff)
            set(MyHandles.ParameterPanel.nDnA.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.nDnA)
            set(MyHandles.ParameterPanel.EAmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EAmax) 
            set(MyHandles.ParameterPanel.EDmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EDmax)
        end
        
        DATASTRUCT.('Constructs').File(i).Data.EA = DATASTRUCT.('Constructs').File(i).Data.YFPfret./DATASTRUCT.('Constructs').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
        DATASTRUCT.('Constructs').File(i).Data.ED = DATASTRUCT.('Constructs').File(i).Data.YFPfret./(DATASTRUCT.('Constructs').File(i).Data.YFPfret+ DATASTRUCT.('Constructs').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);


        
        ND = DATASTRUCT.('Constructs').File(i).Data.CFPdirect./(1-DATASTRUCT.('Constructs').File(i).Data.ED);
        NA = DATASTRUCT.('Constructs').File(i).Data.YFPdirect./(DATASTRUCT.PARAMETERS.gRatio*DATASTRUCT.PARAMETERS.fRatio);        

        
        EA_corrected = DATASTRUCT.('Constructs').File(i).Data.EA - ND*EASlope;
        ED_corrected = DATASTRUCT.('Constructs').File(i).Data.ED - NA*EDSlope; 
        MASK2 = MASK & (ND./NA<DATASTRUCT.MASKPARAMS.NDNA) & (NA./ND<DATASTRUCT.MASKPARAMS.NAND) & (ED_corrected>0);% & (EA_corrected>-0.02)&(ED_corrected>-0.02);   
        EA_corrected_F = EA_corrected(MASK2);
        ED_corrected_F = ED_corrected(MASK2);
        ND_F = ND(MASK2);
        NA_F = NA(MASK2);
        
        if Kd <=0
            Kd = 1000;
        end
        if EAmax <=0 
            EAmax = 0.2;
        end
        if EDmax <=0 
            EDmax = 0.2;
        end
      
        initFitParam = [Kd, EAmax, EDmax];
        
        
        [initFitParam CI] = AutoKDfit(ED_corrected_F,EA_corrected_F,ND_F,NA_F, initFitParam);
        Kd = initFitParam(1);
        EAmax = initFitParam(2);
        EDmax = initFitParam(3);
        
        DATASTRUCT.('Constructs').File(i).FitParams.Kdeff = Kd;
        set(MyHandles.ParameterPanel.Kdeff.Edit,'string',num2str(Kd));
        DATASTRUCT.('Constructs').File(i).FitParams.EAmax = EAmax;
        set(MyHandles.ParameterPanel.EAmax.Edit,'string',num2str(EAmax));
        DATASTRUCT.('Constructs').File(i).FitParams.EDmax = EDmax;
        set(MyHandles.ParameterPanel.EDmax.Edit,'string',num2str(EDmax));
        
        CLIPMATRIX = [initFitParam' CI]; 
        MyFileList = get(MyHandles.Constructs.FileList,'string');
        clipboard('copy',[sprintf('FileName:\t%s\n',MyFileList{i}) sprintf('%s = \t %5.4f \t 95-CI: \t %5.4f \t %5.4f\n','Kd',CLIPMATRIX(1,:),'nD/nA',[1 1 1],'EAmax',CLIPMATRIX(2,:),'EDmax',CLIPMATRIX(3,:))]);
        
    end
    
%     set(gcf,'userdata',DATASTRUCT);  
    UpdatePlots('')
    
    
    
function keypresser(src,evnt)
%     MainFig = gcf;
%     DATASTRUCT = get(src,'userdata');
global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    switch evnt.Key
        case 'c'
            if strcmp(evnt.Modifier{1},'control')
                if strcmp(get(MyHandles.Ax(1),'UserData'),'Constructs')
                    sel = get(MyHandles.Constructs.FileList,'value');
                    Kd = DATASTRUCT.('Constructs').File(sel(1)).FitParams.Kdeff;
                    nDnA = DATASTRUCT.('Constructs').File(sel(1)).FitParams.nDnA; 
                    EAmax = DATASTRUCT.('Constructs').File(sel(1)).FitParams.EAmax;
                    EDmax = DATASTRUCT.('Constructs').File(sel(1)).FitParams.EDmax;
                    clipboard('copy',sprintf('Kd = \t %6.5f \n  nDnA = \t %f \n EAmax = \t %6.5f \n EDmax = \t %6.5f',Kd,nDnA,EAmax,EDmax))
                end
            end
        case 'v'
            if strcmp(evnt.Modifier{1},'control')
                str = clipboard('paste');
                try
                    astr = regexp(str,'\d+\.\d*','match');
                    if length(astr) == 4
                       Kd = str2num(astr{1}); 
                       nDnA = str2num(astr{2}); 
                       EAmax = str2num(astr{3}); 
                       EDmax = str2num(astr{4}); 
                       
                        sel = get(MyHandles.Constructs.FileList,'value');
                        DATASTRUCT.('Constructs').File(sel(1)).FitParams.Kdeff = Kd;
                        DATASTRUCT.('Constructs').File(sel(1)).FitParams.nDnA = nDnA; 
                        DATASTRUCT.('Constructs').File(sel(1)).FitParams.EAmax = EAmax;
                        DATASTRUCT.('Constructs').File(sel(1)).FitParams.EDmax = EDmax;
                        
                        set(MyHandles.ParameterPanel.Kdeff.Edit,'string',Kd)
                        set(MyHandles.ParameterPanel.nDnA.Edit,'string', nDnA)
                        set(MyHandles.ParameterPanel.EAmax.Edit,'string', EAmax) 
                        set(MyHandles.ParameterPanel.EDmax.Edit,'string', EDmax)  
                        
%                         set(MainFig,'userdata',DATASTRUCT)
                        UpdatePlots('')
                    end
                end
            end
            % copy this:
            case 'd'
            if strcmp(evnt.Modifier{1},'control')
                ExportKey('d');
            end
            case 'f'
            if strcmp(evnt.Modifier{1},'control')
                ExportKey('f');
            end
            case 'g'
            if strcmp(evnt.Modifier{1},'control')
                ExportKey('g');
            end
            case 'h'
            if strcmp(evnt.Modifier{1},'control')
                MyFileList = get(MyHandles.Constructs.FileList,'string');
                clipboard('copy',[sprintf('FileName:\t%s\n',MyFileList{get(MyHandles.Constructs.FileList,'value')})]); 
            end
            case 'l'
                AutoLoadFiles();
            
            case 'p'
                autoanalyze();
            case 'a'
                AutoFitCaller();
                % all the way here
    end

function ExportData(varargin)
%     MainFig = gcf;
%     DATASTRUCT = get(MainFig,'userdata');
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    
    SUBDATA = DATASTRUCT.('Constructs').File;
    
    sel = get(MyHandles.Constructs.FileList,'value');
    
    MyColors = jet(length(SUBDATA));
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;
    Gx = PARAMS.gRatio;
    Fx = PARAMS.fRatio;
    DATASTRUCT.PARAMETERS.fRatio = Fx;
    DATASTRUCT.PARAMETERS.gRatio = Gx;
    
    EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
    EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
    if isempty(EASlope) 
        EASlope = 0;
    end
    if isempty(EDSlope) 
        EDSlope = 0;
    end
        
    for i = sel
        CFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.CFP] - PARAMS.BG.CFP;
        YFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.YFP] - PARAMS.BG.YFP;
        FRETTEMP = [DATASTRUCT.('Constructs').File(i).Data.FRET] - PARAMS.BG.FRET;
        
        MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
        MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);
            
        Fin = [CFPTEMP YFPTEMP FRETTEMP]';
        Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';
        
        DATASTRUCT.('Constructs').File(i).Data.CFPdirect  = Sout(:,1);
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect  = Sout(:,2);
        DATASTRUCT.('Constructs').File(i).Data.YFPfret = Sout(:,3);

        DATASTRUCT.('Constructs').File(i).Data.CFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPdirect(MASK<1) = 0;
        DATASTRUCT.('Constructs').File(i).Data.YFPfret(MASK<1) = 0;
    end

    for i = sel
                  
        Kd = DATASTRUCT.('Constructs').File(i).FitParams.Kdeff;% str2num(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));
        nDnA = DATASTRUCT.('Constructs').File(i).FitParams.nDnA;%str2num(get(MyHandles.ParameterPanel.nDnA.Edit,'string'));    
        EAmax = DATASTRUCT.('Constructs').File(i).FitParams.EAmax;%str2num(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));
        EDmax = DATASTRUCT.('Constructs').File(i).FitParams.EDmax; %str2num(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));
        
        if length(sel)==1
            set(MyHandles.ParameterPanel.Kdeff.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.Kdeff)
            set(MyHandles.ParameterPanel.nDnA.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.nDnA)
            set(MyHandles.ParameterPanel.EAmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EAmax) 
            set(MyHandles.ParameterPanel.EDmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EDmax)
        end
        
        DATASTRUCT.('Constructs').File(i).Data.EA = DATASTRUCT.('Constructs').File(i).Data.YFPfret./DATASTRUCT.('Constructs').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
        DATASTRUCT.('Constructs').File(i).Data.ED = DATASTRUCT.('Constructs').File(i).Data.YFPfret./(DATASTRUCT.('Constructs').File(i).Data.YFPfret+ DATASTRUCT.('Constructs').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);


        
        ND = DATASTRUCT.('Constructs').File(i).Data.CFPdirect./(1-DATASTRUCT.('Constructs').File(i).Data.ED);
        NA = DATASTRUCT.('Constructs').File(i).Data.YFPdirect./(DATASTRUCT.PARAMETERS.gRatio*DATASTRUCT.PARAMETERS.fRatio);        
        MASK2 = MASK & (ND./NA<DATASTRUCT.MASKPARAMS.NDNA) & (NA./ND<DATASTRUCT.MASKPARAMS.NAND);
        
        EA_corrected = DATASTRUCT.('Constructs').File(i).Data.EA - ND*EASlope;
        ED_corrected = DATASTRUCT.('Constructs').File(i).Data.ED - NA*EDSlope; 
        
        EA_corrected_F = EA_corrected(MASK2);
        ED_corrected_F = ED_corrected(MASK2);
        ND_F = ND(MASK2);
        NA_F = NA(MASK2);
        
        
        % Compute Dfree / Afree
        Db = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*ND);
        Ab = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*NA);
        
        Dfree = ND.*(1-Db);
        Afree = NA.*(1-Ab);

        

        
        % First Copy Data:
        
        OUTMAT = [(1:sum(MASK2))' Sout(MASK2,:) ND(MASK2,:) NA(MASK2,:) Dfree(MASK2,:) Afree(MASK2,:) EA_corrected(MASK2) ED_corrected(MASK2)]';
        
        DataOutput = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Cell','SCer','SVen','SFRET','ND','NA','Dfree','Afree','EA','ED');
        DataOutputValues = sprintf('%d\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n',OUTMAT);
        clipboard('copy',[DataOutput DataOutputValues])
        
        h = msgbox({'Copied Analysis!' 'Please paste DATA into excel template and, once completed, press OK to continue :)'},'Export FRET Data -- Stage 1', 'help');
        uiwait(h)
        parameternames =  fieldnames(PARAMS);
        parameternames = parameternames(2:end-3);
        Str1 = sprintf('BG.CFP\t%5.4f\nBG.YFP\t%5.4f\nBG.FRET\t%5.4f\nBG.CFP_cutoff\t%5.4f\nBG.YFP_cutoff\t%5.4f\nBG.FRET_cutoff\t%5.4f\n',PARAMS.BG.CFP,PARAMS.BG.YFP,PARAMS.BG.FRET,PARAMS.BG.CFP_CUTOFF,PARAMS.BG.YFP_CUTOFF,PARAMS.BG.FRET_CUTOFF);
        Str2 = sprintf(['%s\\t%%10.9f\\n'],parameternames{:});
        parammatrix = zeros(size(parameternames));
        for kj =1:length(parammatrix)
            parammatrix(kj) = PARAMS.(parameternames{kj});
        end
        clipboard('copy',[Str1 sprintf(Str2,parammatrix)]);
        h = msgbox({'Copied Parameters!' 'Please paste PARAMETERS into excel template and, once completed, press OK to continue :)'},'Export FRET Data -- Stage 1', 'help');
        uiwait(h)
        
        M =  DATASTRUCT.MASKPARAMS;
        Str1 = sprintf('CFP:\t%5.4f\t%5.4f\nYFP:\t%5.4f\t%5.4f\nFRET:\t%5.4f\t%5.4f\n ND/NA:\t %5.4f\n NA/ND\t%5.4f\n',M.CFP.Min,M.CFP.Max,M.YFP.Min,M.YFP.Max,M.FRET.Min,M.FRET.Max,M.NDNA,M.NAND);
        clipboard('copy',Str1);
        h = msgbox({'Copied Mask Parameters!' 'Please paste MASK into excel template and, once completed, press OK to FINISH!!!!'},'Export FRET Data -- Stage 1', 'help');
        uiwait(h)
        
    end
    
    
    %% copy this 
function ExportKey(varargin)
    typeKey = varargin{1};
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;
    
    SUBDATA = DATASTRUCT.('Constructs').File;
    sel = get(MyHandles.Constructs.FileList,'value');
    
    MyColors = jet(length(SUBDATA));
    PARAMS = getPARAMETERS(MyHandles);
    BG  = PARAMS.BG;
    RA1 = PARAMS.RA1;
    RD1 = PARAMS.RD1;
    RD2 = PARAMS.RD2;
    Gx = PARAMS.gRatio;
    Fx = PARAMS.fRatio;
    DATASTRUCT.PARAMETERS.fRatio = Fx;
    DATASTRUCT.PARAMETERS.gRatio = Gx;
    
    EASlope = str2num(get(MyHandles.ParameterPanel.EASlope.Edit,'string'));
    EDSlope = str2num(get(MyHandles.ParameterPanel.EDSlope.Edit,'string'));
    
    
    switch typeKey
        case 'f'
            parameternames =  fieldnames(PARAMS);
            parameternames = parameternames(2:end-3);
            Str1 = sprintf('BG.CFP\t%5.4f\nBG.YFP\t%5.4f\nBG.FRET\t%5.4f\nBG.CFP_cutoff\t%5.4f\nBG.YFP_cutoff\t%5.4f\nBG.FRET_cutoff\t%5.4f\n',PARAMS.BG.CFP,PARAMS.BG.YFP,PARAMS.BG.FRET,PARAMS.BG.CFP_CUTOFF,PARAMS.BG.YFP_CUTOFF,PARAMS.BG.FRET_CUTOFF);
            Str2 = sprintf(['%s\\t%%10.9f\\n'],parameternames{:});
            parammatrix = zeros(size(parameternames));
            for kj =1:length(parammatrix)
                parammatrix(kj) = PARAMS.(parameternames{kj});
            end
            clipboard('copy',[Str1 sprintf(Str2,parammatrix)]);
            return;
        case 'g'
            M =  DATASTRUCT.MASKPARAMS;
            Str1 = sprintf('CFP:\t%5.4f\t%5.4f\nYFP:\t%5.4f\t%5.4f\nFRET:\t%5.4f\t%5.4f\n ND/NA:\t %5.4f\n NA/ND\t%5.4f\n',M.CFP.Min,M.CFP.Max,M.YFP.Min,M.YFP.Max,M.FRET.Min,M.FRET.Max,M.NDNA,M.NAND);
            clipboard('copy',Str1);    
            return;
        case 'd'    
            if isempty(EASlope) 
                EASlope = 0;
            end
            if isempty(EDSlope) 
                EDSlope = 0;
            end

            for i = sel
                CFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.CFP] - PARAMS.BG.CFP;
                YFPTEMP  = [DATASTRUCT.('Constructs').File(i).Data.YFP] - PARAMS.BG.YFP;
                FRETTEMP = [DATASTRUCT.('Constructs').File(i).Data.FRET] - PARAMS.BG.FRET;

                MASK = (CFPTEMP>DATASTRUCT.MASKPARAMS.CFP.Min)&(YFPTEMP>DATASTRUCT.MASKPARAMS.YFP.Min)&(FRETTEMP>DATASTRUCT.MASKPARAMS.FRET.Min);
                MASK = MASK&(CFPTEMP<DATASTRUCT.MASKPARAMS.CFP.Max)&(YFPTEMP<DATASTRUCT.MASKPARAMS.YFP.Max)&(FRETTEMP<DATASTRUCT.MASKPARAMS.FRET.Max);

                Fin = [CFPTEMP YFPTEMP FRETTEMP]';
                Sout = TransformFluorescence(Fin,RA1,RD1,RD2)';

                DATASTRUCT.('Constructs').File(i).Data.CFPdirect  = Sout(:,1);
                DATASTRUCT.('Constructs').File(i).Data.YFPdirect  = Sout(:,2);
                DATASTRUCT.('Constructs').File(i).Data.YFPfret = Sout(:,3);

                DATASTRUCT.('Constructs').File(i).Data.CFPdirect(MASK<1) = 0;
                DATASTRUCT.('Constructs').File(i).Data.YFPdirect(MASK<1) = 0;
                DATASTRUCT.('Constructs').File(i).Data.YFPfret(MASK<1) = 0;
            end

            for i = sel

                Kd = DATASTRUCT.('Constructs').File(i).FitParams.Kdeff;% str2num(get(MyHandles.ParameterPanel.Kdeff.Edit,'string'));
                nDnA = DATASTRUCT.('Constructs').File(i).FitParams.nDnA;%str2num(get(MyHandles.ParameterPanel.nDnA.Edit,'string'));    
                EAmax = DATASTRUCT.('Constructs').File(i).FitParams.EAmax;%str2num(get(MyHandles.ParameterPanel.EAmax.Edit,'string'));
                EDmax = DATASTRUCT.('Constructs').File(i).FitParams.EDmax; %str2num(get(MyHandles.ParameterPanel.EDmax.Edit,'string'));

                if length(sel)==1
                    set(MyHandles.ParameterPanel.Kdeff.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.Kdeff)
                    set(MyHandles.ParameterPanel.nDnA.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.nDnA)
                    set(MyHandles.ParameterPanel.EAmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EAmax) 
                    set(MyHandles.ParameterPanel.EDmax.Edit,'string', DATASTRUCT.('Constructs').File(i).FitParams.EDmax)
                end

                DATASTRUCT.('Constructs').File(i).Data.EA = DATASTRUCT.('Constructs').File(i).Data.YFPfret./DATASTRUCT.('Constructs').File(i).Data.YFPdirect*DATASTRUCT.PARAMETERS.gRatio;
                DATASTRUCT.('Constructs').File(i).Data.ED = DATASTRUCT.('Constructs').File(i).Data.YFPfret./(DATASTRUCT.('Constructs').File(i).Data.YFPfret+ DATASTRUCT.('Constructs').File(i).Data.CFPdirect*DATASTRUCT.PARAMETERS.fRatio);



                ND = DATASTRUCT.('Constructs').File(i).Data.CFPdirect./(1-DATASTRUCT.('Constructs').File(i).Data.ED);
                NA = DATASTRUCT.('Constructs').File(i).Data.YFPdirect./(DATASTRUCT.PARAMETERS.gRatio*DATASTRUCT.PARAMETERS.fRatio);        
                MASK2 = MASK & (ND./NA<DATASTRUCT.MASKPARAMS.NDNA) & (NA./ND<DATASTRUCT.MASKPARAMS.NAND);

                EA_corrected = DATASTRUCT.('Constructs').File(i).Data.EA - ND*EASlope;
                ED_corrected = DATASTRUCT.('Constructs').File(i).Data.ED - NA*EDSlope; 

                EA_corrected_F = EA_corrected(MASK2);
                ED_corrected_F = ED_corrected(MASK2);
                ND_F = ND(MASK2);
                NA_F = NA(MASK2);


                % Compute Dfree / Afree
                Db = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*ND);
                Ab = ((ND+NA+Kd)-sqrt((ND+NA+Kd).^2-4*ND.*NA))./(2*NA);

                Dfree = ND.*(1-Db);
                Afree = NA.*(1-Ab);

                % First Copy Data:

                OUTMAT = [(1:sum(MASK2))' Sout(MASK2,:) ND(MASK2,:) NA(MASK2,:) Dfree(MASK2,:) Afree(MASK2,:) EA_corrected(MASK2) ED_corrected(MASK2)]';

                DataOutput = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Cell','SCer','SVen','SFRET','ND','NA','Dfree','Afree','EA','ED');
                DataOutputValues = sprintf('%d\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n',OUTMAT);
                clipboard('copy',[DataOutput DataOutputValues])
            end
            return;            
    end
    
function AutoLoadFiles(varargin) 
    global DATASTRUCT    
    MyHandles = DATASTRUCT.Handles;
    
    % See if path exists, if not pick a path. 
    if isempty(DATASTRUCT.ActivePath)      
        DATASTRUCT.ActivePath = uigetdir;
        if DATASTRUCT.ActivePath ==0
            errordlg('No Directory Chosen')
            return;
        end
        set(MyHandles.BasePathText,'string',DATASTRUCT.ActivePath);
    end
    
    FFA_Files = dir([DATASTRUCT.ActivePath '/*.ffa']);
    FFA_Names = {FFA_Files.name};
    
    dimer = strncmpi('dimer',FFA_Names,5);
    calib = strncmpi('calib',FFA_Names,5);
    rest = ~(dimer|calib);
    
    % parse calib into background, 
    Calib_FFA_Names = FFA_Names(calib);
    Cer = Calib_FFA_Names(~cellfun(@isempty,strfind(lower(Calib_FFA_Names),'cer')));
    Ven = Calib_FFA_Names(~cellfun(@isempty,strfind(lower(Calib_FFA_Names),'ven')));
    SpFRET = Calib_FFA_Names(~cellfun(@isempty,strfind(lower(Calib_FFA_Names),'sp')));
    
    % blank could be called unstain or blank or bg
    
    BlankIndx = ~cellfun(@isempty,strfind(lower(Calib_FFA_Names),'unst')) | ~cellfun(@isempty,strfind(lower(Calib_FFA_Names),'blank')) | ~cellfun(@isempty,strfind(lower(Calib_FFA_Names),'bg'));
    Blank = Calib_FFA_Names(BlankIndx);
    Dimers = FFA_Names(dimer);
    Constructs = FFA_Names(rest);
    Tags = {'Background','Acceptor','Donor','Dimers','Spurious_FRET','Constructs'};
    
    TableMat = cell(length(FFA_Names),3);
    % First populate Blank
    ibase = 0;
    TableMat((1:length(Blank))+ibase,1) = Blank;
    TableMat((1:length(Blank))+ibase,2) = {'Background'};
    TableMat((1:length(Blank))+ibase,3) = num2cell(ones(length(Blank),1)>0); 
    ibase = ibase + length(Blank);
    
    TableMat((1:length(Ven))+ibase,1) = Ven;
    TableMat((1:length(Ven))+ibase,2) = {'Acceptor'};
    TableMat((1:length(Ven))+ibase,3) = num2cell(ones(length(Ven),1)>0); 
    ibase = ibase + length(Ven);
    
    TableMat((1:length(Cer))+ibase,1) = Cer;
    TableMat((1:length(Cer))+ibase,2) = {'Donor'};
    TableMat((1:length(Cer))+ibase,3) = num2cell(ones(length(Cer),1)>0); 
    ibase = ibase + length(Cer);
    
    TableMat((1:length(Dimers))+ibase,1) = Dimers;
    TableMat((1:length(Dimers))+ibase,2) = {'Dimers'};
    TableMat((1:length(Dimers))+ibase,3) = num2cell(ones(length(Dimers),1)>0); 
    ibase = ibase + length(Dimers);
    
    
    TableMat((1:length(SpFRET))+ibase,1) = SpFRET;
    TableMat((1:length(SpFRET))+ibase,2) = {'Spurious_FRET'};
    TableMat((1:length(SpFRET))+ibase,3) = num2cell(ones(length(SpFRET),1)>0);
    ibase = ibase + length(SpFRET);
    
    TableMat((1:length(Constructs))+ibase,1) = Constructs;
    TableMat((1:length(Constructs))+ibase,2) = {'Constructs'};
    TableMat((1:length(Constructs))+ibase,3) = num2cell(zeros(length(Constructs),1)>0);
    ibase = ibase + length(Constructs);
    
    
    TableMat = TableDlg(TableMat,{'char',Tags,'logical'},[false, true, true],{'FileName','Type','Include?'},{220 'auto' 'auto'});
    if isempty(TableMat)
        return;
    end
    Files2Add = TableMat([TableMat{:,3}]',1:2);
    for jindx = 1:length(Files2Add)
        TempTag = Files2Add{jindx,2};
        TempFile = Files2Add{jindx,1};
        FileList = get(MyHandles.(TempTag).FileList,'string');
        nFiles = length(FileList);
        if (nFiles<1) || (sum(strcmp(FileList,[DATASTRUCT.ActivePath '\' TempFile]))<1)
            FileList{nFiles+1} = [DATASTRUCT.ActivePath '\' TempFile];
            set(MyHandles.(TempTag).FileList,'string',FileList,'value',nFiles+1);
            DATASTRUCT.(TempTag).File(nFiles+1) = loadFFA(FileList{nFiles+1}); 
        end
    end


    function Answer = TableDlg(TableMat,ColFormat,ColEditable,ColName,ColWidth)
    FigAutoload = figure('color',[1 1 1],'units','normalized','UserData','Cancel'); 
    MyTable = uitable('units','normalized','position',[0.1 0.2, 0.8 0.75],'Data',TableMat, 'ColumnFormat',ColFormat,'ColumnEditable',ColEditable,'ColumnName',ColName,'ColumnWidth',ColWidth);
    butOK = uicontrol('units','normalized','position',[0.1 0.05, 0.1, 0.07],'String','OK','callback',@doOK);
    butCancel = uicontrol('units','normalized','position',[0.2 0.05, 0.1, 0.07],'String','Cancel','callback',@doCancel);

    if ishghandle(FigAutoload)
      % Go into uiwait if the figure handle is still valid.
      % This is mostly the case during regular use.
      uiwait(FigAutoload);
    end
    try
    if ishghandle(FigAutoload)
        if strcmp(get(FigAutoload,'UserData'),'OK'),
           Answer = get(MyTable,'Data');
           delete(FigAutoload)
        else
           Answer = {};
           delete(FigAutoload)
        end
    else
        Answer = {};
        delete(FigAutoload)        
    end
    catch
        Answer = {};
    end
    
function doOK(obj, evd)      
    set(gcbf,'UserData','OK');
    uiresume(gcbf);

function doCancel(obj, evd)      
    set(gcbf,'UserData','Cancel');
    uiresume(gcbf);
        
function autoanalyze()
    global DATASTRUCT
    MyHandles = DATASTRUCT.Handles;

    if DATASTRUCT.WORKFLOW.bBGSet == 0
        SetupBackground(1);
    end
    if DATASTRUCT.WORKFLOW.bRASet == 0
        SetupRA(1);
    end
    if DATASTRUCT.WORKFLOW.bRDSet == 0
        SetupRD(1);
    end
    SetupDimers(1);
    
