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
Inputs.Hyd.InitialGeometry = GetInputParameter(C,'GEOMETRY');
Inputs.Hyd.Flow = GetInputParameter(C,'FLOW');
Inputs.Hyd.Slope = GetInputParameter(C,'SLOPE');
Inputs.Hyd.Roughness = GetInputParameter(C,'ROUGHNESS');
Inputs.Hyd.Radius = GetInputParameter(C,'RADIUS');
Inputs.Hyd.DryFlc = GetInputParameter(C,'DRYFLC',0); % drying and flooding threshold [m]
Inputs.Hyd.QTol = GetInputParameter(C,'QTOL',Inputs.Hyd.Flow/1000); % flow tolerance when calculating water level [m3/s]
Inputs.Hyd.Kappa = GetInputParameter(C,'KAPPA',0.4);
Inputs.Hyd.Rho_W = GetInputParameter(C,'RHOW',1000);
Inputs.Hyd.g = GetInputParameter(C,'GRAVITY',9.81);

%% Read in sediment parameters
[Inputs.Sed.SedSize, Inputs.Sed.SedType] = GetInputParameter(C,'SEDSIZE');
Inputs.Sed.Rho_S = GetInputParameter(C,'RHOS',2650);
Inputs.Sed.Porosity = GetInputParameter(C,'POROSITY',0.4);
Inputs.Sed.DA = GetInputParameter(C,'DA'); % Active layer thickness [m]
Inputs.Sed.SedThr = GetInputParameter(C,'SEDTHR',Inputs.Hyd.DryFlc*2); % threshold depth for sediment transport [m]
if Inputs.Sed.SedThr<Inputs.Hyd.DryFlc
    error('SEDTHR must be >= DRYFLC')
end

%% Read in options

% Upwind bedload?
Inputs.Opt.UpwindBedload = GetInputParameter(C,'UPWINDBEDLOAD',1);

% Spiral flow
Inputs.Opt.Spiral.ESpiral = GetInputParameter(C,'ESPIR',1); % Coefficient for effect of spiral flow on bedload transport

% Transport formula
Inputs.Opt.ST.Formula = GetInputParameter(C,'STFORMULA',1);

% Bedslope effects
Inputs.Opt.Slope.Formula = GetInputParameter(C,'BEDSLOPE',0);
switch Inputs.Opt.Slope.Formula
    case 0
        fprintf('No bed slope formulation being used\n')
    case 1
        fprintf('Talmon et al (1995) bed slope calculation\n')
        Inputs.Opt.Slope.A_sh = GetInputParameter(C,'ASH',9);
        Inputs.Opt.Slope.B_sh = GetInputParameter(C,'BSH',0.5);
        Inputs.Opt.Slope.C_sh = GetInputParameter(C,'CSH',0.3);
        Inputs.Opt.Slope.D_sh = GetInputParameter(C,'DSH',0.7);
end

% Bank identification
Inputs.Opt.Bank.ID.Approach = GetInputParameter(C,'BANKID',0);
switch Inputs.Opt.Bank.ID.Approach
    case 0 % no additional inputs required
        fprintf('No bank ID approach being used (i.e. everywhere is a bank)\n')
    case 1
        fprintf('Wet/dry bank identification being used\n')
    case 2
        fprintf('Transporting/non-transporting bank identification being used\n')
    case 3
        fprintf('Bank height bank identification approach being used\n')
        Inputs.Opt.Bank.ID.BHeight = GetInputParameter(C,'BHEIGHT');
    case 4
        fprintf('Bank slope bank identification approach being used\n')
        Inputs.Opt.Bank.ID.BSlope = GetInputParameter(C,'BSLOPE');
end

% Bank stencil
Inputs.Opt.Bank.Stencil.Top = GetInputParameter(C,'BANKTOP',0);
Inputs.Opt.Bank.Stencil.Bottom = GetInputParameter(C,'BANKBOT',0);

% Bank flux calculation
Inputs.Opt.Bank.Flux.Approach = GetInputParameter(C,'BANKFLUX',0);
switch Inputs.Opt.Bank.Flux.Approach
    case 0 % no additional inputs required
        fprintf('No bank erosion flux\n')
    case 1
        fprintf('Excess slope bank erosion flux\n')
        Inputs.Opt.Bank.Flux.Repose = GetInputParameter(C,'REPOSE');
        Inputs.Opt.Bank.Flux.SlipRatio = GetInputParameter(C,'SLIPRATIO',1);
end

% Times
Inputs.Time.dT = GetInputParameter(C,'DT');
Inputs.Time.StartTime = GetInputParameter(C,'STARTTIME',0);
Inputs.Time.EndTime = GetInputParameter(C,'ENDTIME');

% Outputs
Inputs.Outputs.DiagInt = GetInputParameter(C,'DIAGINT',Inputs.Time.dT);
Inputs.Outputs.PlotInt = GetInputParameter(C,'PLOTINT',Inputs.Time.dT);

end

