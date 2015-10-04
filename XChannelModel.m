% Single cross-section 1D (cross-channel) morphological model
% Richard Measures 2015


%% Program setup
addpath('Functions')

%% Read model data and options from intput file
[Inputs] = ReadModelInputs('TestModel.txt','C:\Projects\Research\BankErosionMatlab\Trunk\Inputs\');

%% Prep the model

% Initialise geometry
NCells = size(Inputs.Hyd.InitialGeometry,1);
CellCentrePos = Inputs.Hyd.InitialGeometry(:,1);
CellEdgePos = [CellCentrePos(1) - (CellCentrePos(2) - CellCentrePos(1)) / 2;
               (CellCentrePos(1:end-1) + CellCentrePos(2:end)) / 2;
               CellCentrePos(end) + (CellCentrePos(end) - CellCentrePos(end-1)) / 2];
CellWidth = CellEdgePos(2:end) - CellEdgePos(1:end-1);
BedLevel = Inputs.Hyd.InitialGeometry(:,2);

% Initialise hydraulics
WetCells = ones(NCells,1);

% Initialise sediment
if Inputs.Sed.SedType == 1 % uniform sediment
    Di_m = Inputs.Sed.SedSize;
    Fi = ones(NCells,1);
    BulkFi = Fi;
elseif Inputs.Sed.SedType == 2 % graded sediment
    Di_m = Inputs.Sed.SedSize(:,1)';
    Fi = ones(NCells,1) * Inputs.Sed.SedSize(:,2)';
    BulkFi = ones(NCells,1) * Inputs.Sed.SedSize(:,3)';
end
NFracs = size(Di_m,2);
Di_phi = -log2(Di_m * 1000);
if Inputs.Opt.ST.Formula == 1
    SandFrac = Di_m <= 0.002;
end
SubDg_m = 2.^-sum((ones(NCells,1)*Di_phi) .* BulkFi, 2) ./1000;

% Initialise morpho
Delta_tot = zeros(NCells,1);
Delta_i_tot = zeros(NCells,NFracs);

% Initialise Time
T = Inputs.Time.StartTime - Inputs.Time.dT;

% Set up cross-section plot
%[XsFigureH, BedPlotH, WlPlotH, XsPlotTitleH] = PlotXS(CellCentrePos,CellEdgePos,InitialGeometry(:,2),BedLevel,WL,DryFlc,0);
Dummy = NaN(size(CellCentrePos));
XsFigure = PlotXS(CellCentrePos,CellEdgePos,Inputs.Hyd.InitialGeometry(:,2),BedLevel,NaN,WetCells,Dummy,Dummy,Dummy,Dummy,SubDg_m,0);
PlotT = Inputs.Time.StartTime;
DiagT = Inputs.Time.StartTime;

%% Main Loop

while T < Inputs.Time.EndTime
    
    % Move to next timestep
    T = T + Inputs.Time.dT;
    WetLastTimestep = WetCells;
    
    % Update bed composition
    if Inputs.Sed.SedType == 2
        % NEEDS MORE THOUGHT HERE ABOUT MIXING MODEL...
        % Fi = (Fi * DA + Delta_i_tot * DT) / DA; % Could mix new sed in before moving from active to sub-surface??? Timestep sensitive?
        SubSurfFlux_i = (max(-Delta_tot,0) * ones(1,NFracs)) .* BulkFi - ...
                        (max(Delta_tot,0) * ones(1,NFracs)) .* Fi; % Fractional volumetric flux rate into active layer from subsurface [m3/s/m]
        % Fi = (Fi * DA - SubSurfFlux_i) / DA;
        Fi = (Fi * Inputs.Sed.DA + (Delta_i_tot - SubSurfFlux_i) * Inputs.Time.dT) / Inputs.Sed.DA;
    end

    % Update bed level
    BedLevel = BedLevel + Delta_tot * Inputs.Time.dT / Inputs.Sed.Porosity;
    
    % Calculate basic hydraulics
    H = zeros(NCells,1);
    VAvVel = zeros(NCells,1);
    Tau_S = zeros(NCells,1);
    
    WL = SetQ(CellWidth,BedLevel,WetLastTimestep,Inputs.Hyd);
    WetCells = (WL-BedLevel >= Inputs.Hyd.DryFlc) | (WetLastTimestep & (WL-BedLevel >= Inputs.Hyd.DryFlc / 2));
    H(WetCells) = WL - BedLevel(WetCells);
    VAvVel(WetCells) = H(WetCells).^(2/3) * Inputs.Hyd.Slope^0.5 / Inputs.Hyd.Roughness;
    Tau_S(WetCells) = Inputs.Hyd.Rho_W * Inputs.Hyd.g * VAvVel(WetCells).^2 * Inputs.Hyd.Roughness^2 ./ H(WetCells).^(1/3);

    % Calculate secondary flow
    AlphaSpiral = NaN(NCells,1);
    SpiralIntensity = zeros(NCells,1);
    Tau_N = zeros(NCells,1);
    
    AlphaSpiral(WetCells) = min(sqrt(Inputs.Hyd.g) ./ (Inputs.Hyd.Kappa * (H(WetCells).^(1/6) / Inputs.Hyd.Roughness)), 0.5); % Equation 9.156 in Delft3D-FLOW User Manual or Eq8 Kalkwijk & Booji 1986
    SpiralIntensity(WetCells) = H(WetCells) .* VAvVel(WetCells) / Inputs.Hyd.Radius; % Eq 9.160 (I_be (intensity due to bend) only as I_ce (coriolis) is negligable)
    
    Tau_N(WetCells) = -2 * Inputs.Opt.Spiral.ESpiral * Inputs.Hyd.Rho_W * AlphaSpiral(WetCells).^2 .* (1 - AlphaSpiral(WetCells) / 2) .* VAvVel(WetCells) .* SpiralIntensity(WetCells); % Eq9.145 D3D or Eq26 Kalkwijk & Booji 1986
    Tau_Tot = sqrt(Tau_S.^2 + Tau_N.^2);
    UStar = sqrt(Tau_Tot / Inputs.Hyd.Rho_W);

    % Sediment properties
    Dg_phi = sum((ones(NCells,1)*Di_phi) .* Fi, 2);
    Dg_m = 2.^-Dg_phi ./ 1000;
    SigmaG_phi = sqrt(sum(Fi .* ((ones(NCells,1)*Di_phi) - (Dg_phi*ones(1,NFracs))).^2, 2));

    % Calculate voumetric sediment transport due to flow (in cell centres)
    ActiveCells = (WL - BedLevel) >= Inputs.Sed.SedThr;
    ActiveEdges = [0; ActiveCells(1:end-1) .* ActiveCells(2:end); 0];
    if Inputs.Opt.ST.Formula == 1
        qsiTot_flow = WilcockCrowe(Inputs.Sed.Rho_S, Inputs.Hyd.Rho_W, Inputs.Hyd.g, NCells, NFracs, Di_m, SandFrac, Fi, Dg_m, Tau_Tot, UStar);
    end
    qsiTot_flow(~ActiveCells) = 0;
    qsS_flow_kg = zeros(NCells,1);
    qsS_flow_kg(ActiveCells) = Inputs.Sed.Rho_S * sum(qsiTot_flow(ActiveCells), 2) .* (Tau_S(ActiveCells) ./ Tau_Tot(ActiveCells)); % Streamwise mass transport rate [kg/s/m]
    qsiN_flow = qsiTot_flow .* ((Tau_N ./ Tau_Tot) * ones(1, NFracs));
    qsiN_flow(isnan(qsiN_flow)) = 0;

    % Calculate cell edge parameters
    if Inputs.Opt.UpwindBedload % Upwind bedload
        C2E_Weights = (Tau_N(1:end-1) + Tau_N(2:end)) > 0;
    else % Central
        C2E_Weights = 0.5 * ones(NCells - 1, 1);
    end
    
    qsiTot_flow_edges = Centre2Edge(qsiTot_flow, C2E_Weights);
    qsiTot_flow_edges(~ActiveEdges,:) = 0;
    qsiN_flow_edges = Centre2Edge(qsiN_flow, C2E_Weights);
    qsiN_flow_edges(~ActiveEdges,:) = 0;
    H_edges = Centre2Edge(H, C2E_Weights);
    Tau_Tot_edges = Centre2Edge(Tau_Tot, C2E_Weights);
    Dg_m_edges = Centre2Edge(Dg_m, C2E_Weights);
    
    % Calculate sediment transport due to bed slope
    CellSlope = [0; (BedLevel(2:end) - BedLevel(1:end-1)) ./ (CellCentrePos(2:end) - CellCentrePos(1:end-1)); 0];
    % Note: at the moment, for the purposes of slope analysis we are
    %       assuming unadjusted bedload transport vector is normal to
    %       cross-section and that downstream slope << transverse slope
    %       (i.e. secondary flow neglected, downstream slope neglected for 
    %       bed slope calc). 
    switch Inputs.Opt.Slope.Formula
        case 0
            % No bed slope influence on transport
            qsiN_slope_edges = zeros(NCells+1, NFracs);
        case 1
            % Talmon et al 1995
            ShieldsStressi_edges = (Tau_Tot_edges * ones(1, NFracs)) ./ ((Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * Inputs.Hyd.g * (ones(NCells + 1, 1) * Di_m));
            fTheta = Inputs.Opt.Slope.A_sh * ShieldsStressi_edges.^Inputs.Opt.Slope.B_sh .* ...
                     (((ones(NCells+1, 1) * Di_m)) ./ (H_edges * ones(1, NFracs))).^Inputs.Opt.Slope.C_sh .* ...
                     ((ones(NCells+1,1) * Di_m) ./ (Dg_m_edges * ones(1,NFracs))).^Inputs.Opt.Slope.D_sh;
            qsiN_slope_edges = - qsiTot_flow_edges .* (1 ./ fTheta) .* (CellSlope * ones(1, NFracs));
            qsiN_slope_edges(isnan(qsiN_slope_edges)) = 0;
    end
    
    % Calculate sediment flux rate due to bank erosion
    % Identify banks
    
    % Calculate flux
    
    
    % Erosion/deposition
    Delta_i_flow = qsiN_flow_edges(1:end-1,:) - qsiN_flow_edges(2:end,:); % Fractional volumetric flux rate into cell from neighboring cells due to flow [m3/s/m]
    Delta_i_slope = qsiN_slope_edges(1:end-1,:) - qsiN_slope_edges(2:end,:);
    Delta_i_tot = Delta_i_flow + Delta_i_slope;
    Delta_tot = sum(Delta_i_tot, 2);

    % Outputs
    if T >= PlotT + Inputs.Outputs.PlotInt
        PlotT = PlotT + Inputs.Outputs.PlotInt;
        UpdateXsPlot(XsFigure, BedLevel,WL,WetCells,VAvVel,Tau_S,qsS_flow_kg,Dg_m,SubDg_m,T)
    end
    if T >= DiagT + Inputs.Outputs.DiagInt
        DiagT = DiagT + Inputs.Outputs.DiagInt;
        fprintf('T=%gs, Q=%gm^3/s, q_s=%.2em3/s\n', T, Inputs.Hyd.Flow, sum(qsS_flow_kg))
    end
end

% Final outputs
UpdateXsPlot(XsFigure, BedLevel,WL,WetCells,VAvVel,Tau_S,qsS_flow_kg,Dg_m,SubDg_m,T)
