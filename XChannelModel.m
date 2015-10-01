% Single cross-section 1D (cross-channel) morphological model
% Richard Measures 2015


%% Program setup
addpath('Functions')

%% Read the model input file into a cell array
[FileName,PathName] = uigetfile('*.txt','Select the model input file');
if isequal(FileName,0)
    disp('User selected Cancel')
    return
end
addpath(PathName)
fid = fopen(FileName);
C = textscan(fid, '%[^= ] %*[= ] %s', 'CommentStyle', '%');
fclose(fid);

%% Read in parameters

% Hydraulics
InitialGeometry = GetInputParameter(C,'GEOMETRY');
Flow = GetInputParameter(C,'FLOW');
Slope = GetInputParameter(C,'SLOPE');
Roughness = GetInputParameter(C,'ROUGHNESS');
Radius = GetInputParameter(C,'RADIUS');

% Sediment
[SedSize, SedType] = GetInputParameter(C,'SEDSIZE');
Rho_S = GetInputParameter(C,'RHOS',2650);
Porosity = GetInputParameter(C,'POROSITY',0.4);
STFormula = GetInputParameter(C,'STFORMULA','WC');

% Times
dT = GetInputParameter(C,'DT');
StartTime = GetInputParameter(C,'STARTTIME',0);
EndTime = GetInputParameter(C,'ENDTIME');

% Constants
Rho_W = GetInputParameter(C,'RHOW',1000);
G = GetInputParameter(C,'GRAVITY',9.81);
Kappa = GetInputParameter(C,'KAPPA',0.4);

% Advanced parameters
DryFlc = GetInputParameter(C,'DRYFLC'); % drying and flooding threshold [m]
QTol = GetInputParameter(C,'QTOL',Flow/1000); % flow tolerance when calculating water level [m3/s]

DA = GetInputParameter(C,'DA'); % Active layer thickness [m]
SedThr = GetInputParameter(C,'SEDTHR',DryFlc*2); % threshold depth for sediment transport [m]
if SedThr<DryFlc
    error('SEDTHR must be >= DRYFLC')
end
ESpiral = GetInputParameter(C,'ESPIR',1); % Coefficient for effect of spiral flow on bedload transport

UpwindBedload = GetInputParameter(C,'UPWINDBEDLOAD',1);

% Outputs
DiagInt = GetInputParameter(C,'DIAGINT',dT);
PlotInt = GetInputParameter(C,'PLOTINT',dT);


%% Prep the model

% Initialise geometry
NCells = size(InitialGeometry,1);
CellCentrePos = InitialGeometry(:,1);
CellEdgePos = [CellCentrePos(1) - (CellCentrePos(2) - CellCentrePos(1)) / 2;
               (CellCentrePos(1:end-1) + CellCentrePos(2:end)) / 2;
               CellCentrePos(end) + (CellCentrePos(end) - CellCentrePos(end-1)) / 2];
CellWidth = CellEdgePos(2:end) - CellEdgePos(1:end-1);
BedLevel = InitialGeometry(:,2);

% Initialise hydraulics
WetCells = ones(NCells,1);
WL = NaN;

% Initialise sediment
if SedType == 1 % uniform sediment
    D50i_m = SedSize;
    Fi = ones(NCells,1);
    BulkFi = Fi;
elseif SedType == 2 % graded sediment
    D50i_m = SedSize(:,1)';
    Fi = ones(NCells,1) * SedSize(:,2)';
    BulkFi = ones(NCells,1) * SedSize(:,3)';
end
NFracs = size(D50i_m,2);
D50i_phi = -log2(D50i_m * 1000);
if STFormula == 1
    SandFrac = D50i_m <= 0.002;
end

% Initialise morpho
Delta_tot = zeros(NCells,1);
Delta_i_tot = zeros(NCells,NFracs);

% Initialise Time
T = StartTime - dT;

% Set up cross-section plot
[XsFigureH, BedPlotH, WlPlotH, XsPlotTitleH] = PlotXS(CellCentrePos,CellEdgePos,InitialGeometry(:,2),BedLevel,WL,DryFlc,0);
PlotT = StartTime;
DiagT = StartTime;

%% Main Loop

while T < EndTime
    
    % Move to next timestep
    T = T + dT;
    WetLastTimestep = WetCells;
    
    % Update bed composition
    if SedType == 2
        % NEEDS MORE THOUGHT HERE ABOUT MIXING MODEL...
        % Fi = (Fi * DA + Delta_i_tot * DT) / DA; % Could mix new sed in before moving from active to sub-surface??? Timestep sensitive?
        SubSurfFlux_i = (max(-Delta_tot,0) * ones(1,NFracs)) .* BulkFi - ...
                        (max(Delta_tot,0) * ones(1,NFracs)) .* Fi; % Fractional volumetric flux rate into active layer from subsurface [m3/s/m]
        % Fi = (Fi * DA - SubSurfFlux_i) / DA;
        Fi = (Fi * DA + (Delta_i_tot - SubSurfFlux_i) * dT) / DA;
    end

    % Update bed level
    BedLevel = BedLevel + Delta_tot * dT / Porosity;
    
    % Calculate basic hydraulics
    WL = SetQ(CellWidth,BedLevel,WetLastTimestep,Slope,Roughness,DryFlc,QTol,Flow);
    WetCells = (WL-BedLevel >= DryFlc) | (WetLastTimestep & (WL-BedLevel >= DryFlc / 2));
    H = WL - BedLevel(WetCells);
    VAvVel = H.^(2/3) * Slope^0.5 / Roughness;
    
    Tau_S = zeros(NCells,1);
    Tau_S(WetCells) = Rho_W * G * VAvVel.^2 * Roughness^2 ./ H.^(1/3);

    % Calculate secondary flow
    AlphaSpiral = min(sqrt(G) ./ (Kappa * (H.^(1/6) / Roughness)), 0.5); % Equation 9.156 in Delft3D-FLOW User Manual or Eq8 Kalkwijk & Booji 1986
    SpiralIntensity = H .* VAvVel / Radius; % Eq 9.160 (I_be (intensity due to bend) only as I_ce (coriolis) is negligable)
    
    Tau_N = zeros(NCells,1);
    Tau_N(WetCells) = -2 * ESpiral * Rho_W * AlphaSpiral.^2 .* (1 - AlphaSpiral / 2) .* VAvVel .* SpiralIntensity; % Eq9.145 D3D or Eq26 Kalkwijk & Booji 1986
    Tau_Tot = sqrt(Tau_S.^2 + Tau_N.^2);
    UStar = sqrt(Tau_Tot / Rho_W);

    % Sediment properties
    Dg_phi = sum((ones(NCells,1)*D50i_phi) .* Fi, 2);
    Dg_m = 2.^-Dg_phi / 1000;
    SigmaG_phi = sqrt(sum(Fi .* ((ones(NCells,1)*D50i_phi) - (Dg_phi*ones(1,NFracs))).^2, 2));

    % Calculate sediment transport due to flow
    ActiveCells = (WL - BedLevel) >= SedThr;
    if STFormula == 1
        qsiTot_flow = WilcockCrowe(Rho_S, Rho_W, G, NCells, NFracs, D50i_m, SandFrac, Fi, Dg_m, Tau_Tot, UStar);
    end
    qsiTot_flow(~ActiveCells) = 0;
    qsS = sum(sum(qsiTot_flow(ActiveCells), 2) .* (Tau_S(ActiveCells) ./ Tau_Tot(ActiveCells))); % Total volumetric transport rate through section
    qsiN_flow = qsiTot_flow .* ((Tau_N ./ Tau_Tot) * ones(1, NFracs));
    qsiN_flow(isnan(qsiN_flow)) = 0;

    % Calculate cell edge sediment flux rates
    if UpwindBedload % Upwind bedload
        C2E_Weights = (Tau_N(1:end-1) + Tau_N(2:end)) > 0;
    else % Central
        C2E_Weights = 0.5 * ones(NCells - 1, 1);
    end

    qsiN_flow_edges = [zeros(1,NFracs) ;
                       qsiN_flow(1:end-1,:) .* (C2E_Weights*ones(1, NFracs)) + qsiN_flow(2:end,:) .* (1-(C2E_Weights*ones(1, NFracs))); 
                       zeros(1,NFracs)];

    % Erosion/deposition
    Delta_i_flow = qsiN_flow_edges(1:end-1,:) - qsiN_flow_edges(2:end,:); % Fractional volumetric flux rate into cell from neighboring cells due to flow [m3/s/m]
    Delta_i_tot = Delta_i_flow;
    Delta_tot = sum(Delta_i_tot, 2);

    % Outputs
    if T >= PlotT + PlotInt
        PlotT = PlotT + PlotInt;
        UpdateXsPlot(BedPlotH, WlPlotH, XsPlotTitleH, BedLevel,WL,DryFlc,T)
    end
    if T >= DiagT + DiagInt
        DiagT = DiagT + DiagInt;
        fprintf('T=%gs, Q=%gm^3/s, q_s=%.2em3/s\n', T, Flow, qsS)
    end
end

% Final outputs
UpdateXsPlot(BedPlotH, WlPlotH, XsPlotTitleH, BedLevel,WL,DryFlc,T)
