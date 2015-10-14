function [Inputs] = ReadModelInputs(FileName,PathName)
% Read in XChannelModel model input file to structure array
% [Inputs] = ReadModelInputs(FileName,PathName)
% both inputs are optional

%% Get file name and path if not specified as input
if ~exist('FileName','var')
    [FileName,PathName] = uigetfile('*.txt','Select the model input file');
    if isequal(FileName,0)
        error('User selected Cancel')
    end
end

if exist('PathName','var')
    addpath(PathName)
end

%% Read the model input file into a cell array
fid = fopen(FileName);
C = textscan(fid, '%[^= ] %*[= ] %s', 'CommentStyle', '%');
fclose(fid);

%% Read in hydraulics parameters
Inputs.Hyd.InitialGeometry = GetInputParameter(C,'Geometry');
[Inputs.Hyd.Flow, Inputs.Hyd.FlowType] = GetInputParameter(C,'Flow');
if Inputs.Hyd.FlowType == 2 % Flow from timeseries file
    Inputs.Hyd.FlowTS = Inputs.Hyd.Flow;
    Inputs.Hyd.Flow = Inputs.Hyd.FlowTS(1,2);
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
[Inputs.Sed.SedSize, Inputs.Sed.SedType] = GetInputParameter(C,'SedSize');
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
switch Inputs.Opt.Bank.ID.Approach
    case 0 % no additional inputs required
        fprintf('No bank ID approach being used (i.e. everywhere is a bank)\n')
    case 1
        fprintf('Wet/dry bank identification being used\n')
    case 2
        fprintf('Transporting/non-transporting bank identification being used\n')
    case 3
        fprintf('Bank height bank identification approach being used\n')
        Inputs.Opt.Bank.ID.BHeight = GetInputParameter(C,'BHeight');
    case 4
        fprintf('Bank slope bank identification approach being used\n')
        Inputs.Opt.Bank.ID.BSlope = GetInputParameter(C,'BSlope');
end

% Bank stencil
Inputs.Opt.Bank.Stencil.Top = GetInputParameter(C,'BankTop',0);
switch Inputs.Opt.Bank.Stencil.Top
    case 0
        fprintf('No bank top stencil being used - adjacent cells\n')
end
if Inputs.Opt.Bank.Stencil.Top ~=0
    Inputs.Opt.Bank.Stencil.TopCellLim = GetInputParameter(C,'TopCellLim',1);
    Inputs.Opt.Bank.Stencil.TopDistLim = GetInputParameter(C,'TopDistLim',1);
end

Inputs.Opt.Bank.Stencil.Bottom = GetInputParameter(C,'BankBot',0);
switch Inputs.Opt.Bank.Stencil.Bottom
    case 0
        fprintf('No bank bottom stencil being used - adjacent cells\n')
    case 1
        fprintf('Maximum slope curvature bottom stencil being used - adjacent cells\n')
end
if Inputs.Opt.Bank.Stencil.Bottom ~=0
    Inputs.Opt.Bank.Stencil.BotCellLim = GetInputParameter(C,'BotCellLim',1);
    Inputs.Opt.Bank.Stencil.BotDistLim = GetInputParameter(C,'BotDistLim',1);
end

% Bank flux calculation
Inputs.Opt.Bank.Flux.Approach = GetInputParameter(C,'BankFlux',0);
switch Inputs.Opt.Bank.Flux.Approach
    case 0 % no additional inputs required
        fprintf('No bank erosion flux\n')
    case 1
        fprintf('Excess slope bank erosion flux\n')
        Inputs.Opt.Bank.Flux.Repose = GetInputParameter(C,'Repose');
        Inputs.Opt.Bank.Flux.SlipRatio = GetInputParameter(C,'SlipRatio',1);
    case 2
        fprintf('Bank erosion flux proportional to bank toe erosion\n')
        Inputs.Opt.Bank.Flux.ThetSD = GetInputParameter(C,'ThetSD',0.5);
    case 3
        fprintf('Bank erosion flux proportional to bank toe transport rate\n')
        Inputs.Opt.Bank.Flux.QsBeRatio = GetInputParameter(C,'QsBeRatio');
    case 4
        fprintf('Bank erosion flux proportional to (bank toe transport rate * slope)\n')
        Inputs.Opt.Bank.Flux.BErodibility = GetInputParameter(C,'BErodibility');
end
Inputs.Opt.Bank.Flux.StencilMix = GetInputParameter(C,'StencilMix',0);
Inputs.Opt.Bank.Flux.StoredBE = GetInputParameter(C,'StoredBE',0);

% Other bank erosion options
Inputs.Opt.Bank.Flux.Approach = GetInputParameter(C,'BankFlux',0);

% Times
Inputs.Time.dT = GetInputParameter(C,'dT');
Inputs.Time.StartTime = GetInputParameter(C,'StartTime',0);
Inputs.Time.EndTime = GetInputParameter(C,'EndTime');

% Outputs
Inputs.Outputs.DiagInt = GetInputParameter(C,'DiagInt',Inputs.Time.dT);
Inputs.Outputs.PlotInt = GetInputParameter(C,'PlotInt',Inputs.Time.dT);

end

