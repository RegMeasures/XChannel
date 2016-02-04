function [Inputs] = ReadModelInputs(FileName)
% Read in XChannelModel model input file to structure array
% [Inputs] = ReadModelInputs(FileName,PathName)
% 

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
Inputs.Hyd.ManningN = GetInputParameter(C,'ManningN');
Inputs.Hyd.Radius = GetInputParameter(C,'Radius');

Inputs.Hyd.ESpiral = GetInputParameter(C,'ESpir',1);                       % Coefficient for effect of spiral flow on bedload transport
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

%% Bedload transport due to flow
Inputs.ST.UpwindBedload = GetInputParameter(C,'UpwindBedload',1);         % Scheme for calculating bedload at cell edges: 1 = upwind, 2 = central
Inputs.ST.Formula = GetInputParameter(C,'STFormula',1);                % Sediment transport formula
Inputs.ST.ThetaCrit = GetInputParameter(C,'ThetaCrit',0.047);
Inputs.ST.MPMcoef = GetInputParameter(C,'MPMcoef',4.93);
Inputs.ST.MPMexponent = GetInputParameter(C,'MPMexponent',1.6);
switch Inputs.ST.Formula
    case 0
        fprintf('No streamwise bedload transport\n')
    case 1
        fprintf('Meyer-Peter-Muller bedload transport formula\n')
        fprintf('MPM bedload formula coefficients: ThetaCrit = %g, MPMcoef = %g, MPMexponent = %g\n',...
                Inputs.ST.ThetaCrit, Inputs.ST.MPMcoef,...
                Inputs.ST.MPMexponent)
    case 2
        fprintf('Wilcock-Crowe bedload transport formula\n')
end
% Hiding and exposure correction for graded sediment
Inputs.ST.HidExp = GetInputParameter(C,'HidExp',0);
Inputs.ST.Gamma = GetInputParameter(C,'Gamma');
switch Inputs.ST.HidExp
    case 0
        fprintf('No hiding function being used\n')
    case 1
        fprintf('Ashida and Michiue hiding function selected\n')
    case 2
        fprintf('Wilcock and Crowe hiding function selected\n')
    case 3
        fprintf('Parker Klingeman and McLean hiding function selected\n')
        fprintf('Hiding function exponent (Gamma) = %g\n', Inputs.ST.Gamma)
end

%% Bedslope effects on transport
Inputs.Slope.Formula = GetInputParameter(C,'BedSlope',0);
Inputs.Slope.BetaStar = GetInputParameter(C,'BetaStar');
Inputs.Slope.m = GetInputParameter(C,'m');
Inputs.Slope.A_sh = GetInputParameter(C,'Ash',9);
Inputs.Slope.B_sh = GetInputParameter(C,'Bsh',0.5);
Inputs.Slope.C_sh = GetInputParameter(C,'Csh',0.3);
Inputs.Slope.D_sh = GetInputParameter(C,'Dsh',0.7);
switch Inputs.Slope.Formula
    case 0
        fprintf('No bed slope formulation being used\n')
    case 1
        fprintf('General bed slope formulation being used (Sekine and Parker)\n')
        fprintf('Bed slope formulation coefficients: BetaStar = %g, m = %g\n',...
                Inputs.Slope.BetaStar, Inputs.Slope.m);
    case 2
        fprintf('Talmon et al (1995) bed slope calculation\n')
        fprintf('Talmon et al coefficients: Ash = %g, Bsh = %g, Csh = %g, Dsh = %g',...
                Inputs.Slope.A_sh, Inputs.Slope.B_sh,...
                Inputs.Slope.C_sh, Inputs.Slope.D_sh);
end

%% Bank erosion
% Bank identification
Inputs.Bank.ID.Approach = GetInputParameter(C,'BankID',0);
Inputs.Bank.ID.BHeight = GetInputParameter(C,'BHeight');
Inputs.Bank.ID.BSlope = GetInputParameter(C,'BSlope');
switch Inputs.Bank.ID.Approach
    case 0 % no additional inputs required
        fprintf('No bank ID approach being used (i.e. everywhere is a bank)\n')
    case 1
        fprintf('Wet/dry bank identification\n')
    case 2
        fprintf('Transporting/non-transporting bank identification\n')
    case 3
        fprintf('Bank height bank identification approach\n')
        fprintf('Threshold bank height for bank identification (BHeight) = %g\n',...
                Inputs.Bank.ID.BHeight)
    case 4
        fprintf('Bank slope bank identification approach\n')
        fprintf('Threshold bank slope for bank identification (BSlope) = %g\n',...
                Inputs.Bank.ID.BHeight)
end

% Bank stencil
Inputs.Bank.Stencil.Top = GetInputParameter(C,'BankTop',0);
Inputs.Bank.Stencil.TopCellLim = GetInputParameter(C,'TopCellLim',1);
Inputs.Bank.Stencil.TopDistLim = GetInputParameter(C,'TopDistLim',9999);
switch Inputs.Bank.Stencil.Top
    case 0
        fprintf('No bank top stencil being used - adjacent cells\n')
end

Inputs.Bank.Stencil.Bottom = GetInputParameter(C,'BankBot',0);
Inputs.Bank.Stencil.BotCellLim = GetInputParameter(C,'BotCellLim',1);
Inputs.Bank.Stencil.BotDistLim = GetInputParameter(C,'BotDistLim',9999);
switch Inputs.Bank.Stencil.Bottom
    case 0
        fprintf('No bank bottom stencil being used - adjacent cells\n')
    case 1
        fprintf('Maximum slope curvature bottom stencil being used\n')
    case 2
        fprintf('Lowest elevation cell bottom stencil being used\n')
end

% Bank trigger
Inputs.Bank.Trigger.BTrigger = GetInputParameter(C,'BTrigger',0);
Inputs.Bank.Trigger.BTHeight = GetInputParameter(C,'BTHeight');
Inputs.Bank.Trigger.BTSlope = GetInputParameter(C,'BTSlope');
switch Inputs.Bank.Trigger.BTrigger
    case 0
        fprintf('No bank trigger approach (i.e. every bank is active)\n')
    case 1
        fprintf('Threshold height bank trigger selected\n')
        fprintf('Threshold height for bank trigger (BTHeight) = %g\n',...
                Inputs.Bank.Trigger.BTHeight)
    case 2
        fprintf('Degrading toe bank trigger selected\n')
    case 3
        fprintf('Threshold slope bank trigger selected\n')
        fprintf('Threshold slope for bank trigger (BTSlope) = %g\n',...
                Inputs.Bank.Trigger.BTSlope)
end

% Bank flux calculation
Inputs.Bank.Flux.Approach = GetInputParameter(C,'BankFlux',0);
Inputs.Bank.Flux.Repose = GetInputParameter(C,'Repose');
Inputs.Bank.Flux.SlipRatio = GetInputParameter(C,'SlipRatio',1);
Inputs.Bank.Flux.ThetSD = GetInputParameter(C,'ThetSD',0.5);
Inputs.Bank.Flux.QsBeRatio = GetInputParameter(C,'QsBeRatio');
Inputs.Bank.Flux.BErodibility = GetInputParameter(C,'BErodibility');
switch Inputs.Bank.Flux.Approach
    case 0 % no additional inputs required
        fprintf('No bank erosion flux\n')
    case 1
        fprintf('Excess slope bank erosion flux\n')
        fprintf('Limiting slope (Repose) = %g, Excess slipped per timestep (SlipRatio) = %g\n',...
                Inputs.Bank.Flux.Repose, Inputs.Bank.Flux.SlipRatio)
    case 2
        fprintf('Bank erosion flux proportional to bank toe erosion\n')
        fprintf('Proportion of toe erosion redistributed to bank top (ThetSD) = %g\n',...
                Inputs.Bank.Flux.ThetSD)
    case 3
        fprintf('Bank erosion flux proportional to bank toe transport rate\n')
        fprintf('Ratio of bank erosion flux to bank toe transport rate (QsBeRatio) = %g\n',...
                Inputs.Bank.Flux.QsBeRatio)
    case 4
        fprintf('Bank erosion flux proportional to (bank toe transport rate * slope)\n')
        fprintf('Ratio of bank erosion flux to bank toe transport rate * slope (BErodibility) = %g\n',...
                Inputs.Bank.Flux.BErodibility)
end

Inputs.Bank.Flux.StencilMix = GetInputParameter(C,'StencilMix',0);
if (Inputs.Bank.Stencil.Top > 0 || Inputs.Bank.Stencil.Bottom > 0) && Inputs.Sed.SedType == 2
    switch Inputs.Bank.Flux.StencilMix
        case 0
            fprintf('No mixing across stencil (i.e. material is transported from bank top directly to bank toe)')
        case 1
            fprintf('Mix sediment through active layer when transporting across stencil')
    end
end

% Bank updating
Inputs.Bank.Update.StoredBE = GetInputParameter(C,'StoredBE',0);
switch Inputs.Bank.Update.StoredBE
    case 0
        fprintf('Continuous bank erosion')
    case 1
        fprintf('Bank erosion stored until bank top cell can be dropped to toe elevation')
end

%% Times
Inputs.Time.dT = GetInputParameter(C,'dT');
Inputs.Time.StartTime = GetInputParameter(C,'StartTime',0);
Inputs.Time.EndTime = GetInputParameter(C,'EndTime');

%% Outputs
Inputs.Outputs.DiagInt = GetInputParameter(C,'DiagInt',Inputs.Time.dT);
Inputs.Outputs.PlotInt = GetInputParameter(C,'PlotInt',Inputs.Time.dT);
Inputs.Outputs.VideoOut = GetInputParameter(C,'VideoOut',0);

end

