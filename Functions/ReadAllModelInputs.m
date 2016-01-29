function [Inputs] = ReadAllModelInputs(FileName)
% Read in XChannelModel model input file to structure array
% [Inputs] = ReadModelInputs(FileName,PathName)
% 
% Same as ReadModelInputs except it reads in many optional parameters
% irrespective of whether they are required. This is to enable changing of
% options in SensitivityAnalysis

%% Get file name and path if not specified as input
if ~exist('FileName','var')
    [FileName,FilePath] = uigetfile('*.txt','Select the model input file');
    if isequal(FileName,0)
        error('User selected Cancel')
    end
    FileName = fullfile(FilePath,FileName);
end

%% Read the model input file into a cell array

Inputs.FileName = FileName;
FilePath = fileparts(FileName);

fid = fopen(FileName);
C = textscan(fid, '%[^= ] %*[= ] %s', 'CommentStyle', '%');
fclose(fid);

%% Read in hydraulics parameters

% Geometry
[Inputs.Hyd.InitialGeometry, Type] = GetInputParameter(C,'Geometry',[],FilePath);
if Type ~= 2
    error('Geometry file %s not found',Inputs.Hyd.InitialGeometry)
end

% Flow
[Inputs.Hyd.Flow, Inputs.Hyd.FlowType] = GetInputParameter(C,'Flow',[],FilePath);
if Inputs.Hyd.FlowType == 2 % Flow from timeseries file
    Inputs.Hyd.FlowTS = Inputs.Hyd.Flow;
    Inputs.Hyd.Flow = Inputs.Hyd.FlowTS(1,2);
elseif Inputs.Hyd.FlowType == 3 % Missing file or badly formed number?
    error('Flow file %s not found',Inputs.Hyd.Flow)
end

Inputs.Hyd.Slope = GetInputParameter(C,'Slope');
Inputs.Hyd.Roughness = GetInputParameter(C,'Roughness');
Inputs.Hyd.Radius = GetInputParameter(C,'Radius');
Inputs.Hyd.DryFlc = GetInputParameter(C,'DryFlc',0);                       % drying and flooding threshold [m]
Inputs.Hyd.QTol = GetInputParameter(C,'QTol',0.001);                       % flow tolerance when calculating water level [m3/s]
Inputs.Hyd.ItMax = GetInputParameter(C,'ItMax',20);                        % maximum number of iterations for setting WL
Inputs.Hyd.Kappa = GetInputParameter(C,'Kappa',0.4);
Inputs.Hyd.Rho_W = GetInputParameter(C,'RhoW',1000);
Inputs.Hyd.g = GetInputParameter(C,'Gravity',9.81);

%% Read in sediment parameters
[Inputs.Sed.SedSize, Inputs.Sed.SedType] = GetInputParameter(C,'SedSize',[],FilePath);
if Inputs.Sed.SedType == 3;
    error('Sediment file %s not found',Inputs.Sed.SedSize)
end

Inputs.Sed.Rho_S = GetInputParameter(C,'RhoS',2650);
Inputs.Sed.Porosity = GetInputParameter(C,'Porosity',0.4);
Inputs.Sed.DA = GetInputParameter(C,'DA');                                 % Active layer thickness [m]
Inputs.Sed.SedThr = GetInputParameter(C,'SedThr',Inputs.Hyd.DryFlc*2);     % threshold depth for sediment transport [m]
if Inputs.Sed.SedThr<Inputs.Hyd.DryFlc
    error('SedThr must be >= DryFlc')
end

%% Read in options

% Upwind bedload?
Inputs.Opt.UpwindBedload = GetInputParameter(C,'UpwindBedload',1);

% Spiral flow
Inputs.Opt.Spiral.ESpiral = GetInputParameter(C,'ESpir',1); % Coefficient for effect of spiral flow on bedload transport

% Transport formula
Inputs.Opt.ST.Formula = GetInputParameter(C,'STFormula',1);
switch Inputs.Opt.ST.Formula
    case 0
        fprintf('No bedload transport formula being used\n')
    case 1
        fprintf('Meyer-Peter-Muller bedload transport formula being used\n')
        Inputs.Opt.ST.ThetaCrit = GetInputParameter(C,'ThetaCrit',0.047);
        Inputs.Opt.ST.MPMcoef = GetInputParameter(C,'MPMcoef',4.93);
        Inputs.Opt.ST.MPMexponent = GetInputParameter(C,'MPMexponent',1.6);
    case 2
        fprintf('Wilcock-Crowe bedload transport formula being used\n')
end

% Hiding and exposure correction
Inputs.Opt.ST.HidExp = GetInputParameter(C,'HidExp',0);
switch Inputs.Opt.ST.HidExp
    case 0
        fprintf('No hiding function being used\n')
    case 1
        fprintf('Ashida and Michiue hiding function selected\n')
    case 2
        fprintf('Wilcock and Crowe hiding function selected\n')
    case 3
        fprintf('Parker Klingeman and McLean hiding function selected\n')
        Inputs.Opt.ST.Gamma = GetInputParameter(C,'Gamma');
end

% Bedslope effects
Inputs.Opt.Slope.Formula = GetInputParameter(C,'BedSlope',0);
switch Inputs.Opt.Slope.Formula
    case 0
        fprintf('No bed slope formulation being used\n')
    case 1
        fprintf('General bed slope formulation being used (Sekine and Parker)\n')
        Inputs.Opt.Slope.BetaStar = GetInputParameter(C,'BetaStar');
        Inputs.Opt.Slope.m = GetInputParameter(C,'m');
    case 2
        fprintf('Talmon et al (1995) bed slope calculation\n')
        Inputs.Opt.Slope.A_sh = GetInputParameter(C,'Ash',9);
        Inputs.Opt.Slope.B_sh = GetInputParameter(C,'Bsh',0.5);
        Inputs.Opt.Slope.C_sh = GetInputParameter(C,'Csh',0.3);
        Inputs.Opt.Slope.D_sh = GetInputParameter(C,'Dsh',0.7);
end

% Bank identification
Inputs.Opt.Bank.ID.Approach = GetInputParameter(C,'BankID',0);
Inputs.Opt.Bank.ID.BHeight = GetInputParameter(C,'BHeight');
Inputs.Opt.Bank.ID.BSlope = GetInputParameter(C,'BSlope');

% Bank stencil
Inputs.Opt.Bank.Stencil.Top = GetInputParameter(C,'BankTop',0);
Inputs.Opt.Bank.Stencil.TopCellLim = GetInputParameter(C,'TopCellLim',1);
Inputs.Opt.Bank.Stencil.TopDistLim = GetInputParameter(C,'TopDistLim',9999);

Inputs.Opt.Bank.Stencil.Bottom = GetInputParameter(C,'BankBot',0);
Inputs.Opt.Bank.Stencil.BotCellLim = GetInputParameter(C,'BotCellLim',1);
Inputs.Opt.Bank.Stencil.BotDistLim = GetInputParameter(C,'BotDistLim',9999);

% Bank trigger
Inputs.Opt.Bank.Trigger.BTrigger = GetInputParameter(C,'BTrigger',0);
Inputs.Opt.Bank.Trigger.BTHeight = GetInputParameter(C,'BTHeight');
Inputs.Opt.Bank.Trigger.BTSlope = GetInputParameter(C,'BTSlope');

% Bank flux calculation
Inputs.Opt.Bank.Flux.Approach = GetInputParameter(C,'BankFlux',0);
Inputs.Opt.Bank.Flux.Repose = GetInputParameter(C,'Repose');
Inputs.Opt.Bank.Flux.SlipRatio = GetInputParameter(C,'SlipRatio',1);
Inputs.Opt.Bank.Flux.ThetSD = GetInputParameter(C,'ThetSD',0.5);
Inputs.Opt.Bank.Flux.QsBeRatio = GetInputParameter(C,'QsBeRatio');
Inputs.Opt.Bank.Flux.BErodibility = GetInputParameter(C,'BErodibility');
Inputs.Opt.Bank.Flux.StencilMix = GetInputParameter(C,'StencilMix',0);

% Bank updating
Inputs.Opt.Bank.Update.StoredBE = GetInputParameter(C,'StoredBE',0);

%% Times
Inputs.Time.dT = GetInputParameter(C,'dT');
Inputs.Time.StartTime = GetInputParameter(C,'StartTime',0);
Inputs.Time.EndTime = GetInputParameter(C,'EndTime');

%% Outputs
Inputs.Outputs.DiagInt = GetInputParameter(C,'DiagInt',Inputs.Time.dT);
Inputs.Outputs.PlotInt = GetInputParameter(C,'PlotInt',Inputs.Time.dT);
Inputs.Outputs.VideoOut = GetInputParameter(C,'VideoOut',0);

end

