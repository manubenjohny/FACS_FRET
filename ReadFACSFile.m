function [FACSData] = ReadFACSFile(fname)
%%% Function reads FCS files for FRET 
try 
    [NumericArray, TextArray] = fca_readfcs(fname);
catch
    disp('Error Reading File.')
    return;    
end

HEADER.FileName = TextArray.filename;
HEADER.FilePath = TextArray.filepath;
HEADER.cytometry = TextArray.cytometry;
HEADER.system = TextArray.system;


NewHeader.FileName = HEADER.FileName;
NewHeader.FilePath = HEADER.FilePath;
NewHeader.cytometry = HEADER.cytometry;
NewHeader.system = HEADER.system;

if TextArray.NumOfPar == size(NumericArray,2)      
   Names = regexprep({TextArray.par.name},'[-:,/#$&*><|~@! ]', ''); 
   for  i = 1: size(NumericArray, 2)
        TempFACS.(Names{i}) = NumericArray(:,i);
        HEADER.(Names{i}).PMTV = TextArray.par(i).PMTV;
        
   end   
end


switch HEADER.cytometry
    case 'Attune Cytometric Software v2.1.0.8626 (BV)'
        CalibrationFile = './PMTcalib.mat';
        % Attune
        MyFields     = {'FSCA','FSCH','SSCA','SSCH','VL1A','BL1A','VL2A'};
        MyOutputName = {'FSCA','FSCH','SSCA','SSCH','CFP' ,'YFP' ,'FRET'};
        PMTcalibFields = {'FSCHc','FSCHc','SSCHc','SSCHc','VL1Ac' ,'BL1Ac' ,'VL2Ac'};
        
        STDGAIN      = {  1750,  1750,  2250,  2250,  1000,  1000,  1000};
        PMTcalib = load(CalibrationFile);
        for i = 1:length(MyFields)
%             FCSData.Header.(MyOutputName{i}).PMTV = Temp.Header.(MyFields{i}).PMTV;
            defaultPMTV = STDGAIN{i};
            currentPMTV = str2num(HEADER.(MyFields{i}).PMTV);
            xi = [currentPMTV defaultPMTV];
            YCalib = PMTcalib.(PMTcalibFields{i});
            XCalib = PMTcalib.PMTVc;
            NANFilter = (isnan(YCalib)|isnan(XCalib));
            YCalib = YCalib(~NANFilter);
            XCalib = XCalib(~NANFilter);
            yi = exp(interp1(log(XCalib),log(YCalib),log(xi)));
            scalefactor = yi(2)/yi(1);
            OUTPUT.(MyOutputName{i}) = TempFACS.(MyFields{i})*scalefactor;    
            NewHeader.(MyOutputName{i}).PMTV = defaultPMTV;
            NewHeader.(MyOutputName{i}).AcqDye = MyFields{i};
            NewHeader.(MyOutputName{i}).AcqPMTV = currentPMTV;
        end
        
    case 'LSRFortessa' 
        % Attune
        MyFields     = {'FSCA','FSCH','SSCA','SSCH','PacificBlueA','FITCA','AmCyanA'};
        MyOutputName = {'FSCA','FSCH','SSCA','SSCH','CFP' ,'YFP' ,'FRET'};
        STDGAIN      = {210,  210,  240,  240,  270,  210,  230};
        Gs = [100:25:400];
        FITC = [2.682, 15.8363, 45.3024,108.8512, 285.8597, 736.9433, 1758.0077,3437.6977, 6530.2483, 11383.235,18925.869,30617.621,50571.2604];
        for i = 1:length(MyFields)
            defaultPMTV = STDGAIN{i};
            currentPMTV = str2num(HEADER.(MyFields{i}).PMTV);
            xi = [currentPMTV defaultPMTV];
            yi = exp(interp1(log(Gs),log(FITC),log(xi)));
            scalefactor = yi(2)/yi(1);
            OUTPUT.(MyOutputName{i}) = TempFACS.(MyFields{i})*scalefactor;  
            NewHeader.(MyOutputName{i}).PMTV = defaultPMTV;
            NewHeader.(MyOutputName{i}).AcqDye = MyFields{i};
            NewHeader.(MyOutputName{i}).AcqPMTV = currentPMTV;
        end         
    case 'LSRII'
               % Attune
        MyFields     = {'FSCA','FSCH','SSCA','SSCH','BV421A','FITCA','BV510A'};
        MyOutputName = {'FSCA','FSCH','SSCA','SSCH','CFP' ,'YFP' ,'FRET'};
        PMTcalibFields = {'FSCHc','FSCHc','SSCHc','SSCHc','PacificBlueAc' ,'FITCAc' ,'AmCyanAc'};
        STDGAIN      = {416,  416,  302,  302,  300,  300,  300};
        Gs = [100:25:400];
        FITC = [2.682, 15.8363, 45.3024,108.8512, 285.8597, 736.9433, 1758.0077,3437.6977, 6530.2483, 11383.235,18925.869,30617.621,50571.2604];
        for i = 1:length(MyFields)
            defaultPMTV = STDGAIN{i};
            currentPMTV = str2num(HEADER.(MyFields{i}).PMTV);
            xi = [currentPMTV defaultPMTV];
            yi = exp(interp1(log(Gs),log(FITC),log(xi)));
            scalefactor = 1;%yi(2)/yi(1)
            OUTPUT.(MyOutputName{i}) = TempFACS.(MyFields{i})*scalefactor;  
            NewHeader.(MyOutputName{i}).PMTV = defaultPMTV;
            NewHeader.(MyOutputName{i}).AcqDye = MyFields{i};
            NewHeader.(MyOutputName{i}).AcqPMTV = currentPMTV;
        end   
        
        
        
        
        
end

FACS = OUTPUT;

FACSData.Data = FACS;
FACSData.Header = NewHeader;

% HEADER.Path = Tex


%  Temp = load(MyFile,'-mat','Header','FILTData');
