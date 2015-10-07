% Single cross-section 1D (cross-channel) morphological model
% Richard Measures 2015


%% Program setup
addpath('Functions')

%% Read model data and options from intput file
[Inputs] = ReadModelInputs('TestModel.txt','C:\Projects\Research\BankErosionMatlab\Trunk\Inputs\');

%% Initialise the main variable structs
[Cell, Edge, Frac, Bank] = InitialiseVariables(Inputs);

%% Set up cross-section plot
XsFigure = PlotXS(Cell,Edge,Bank,Inputs.Hyd.InitialGeometry(:,2),NaN,0);
PlotT = Inputs.Time.StartTime;
DiagT = Inputs.Time.StartTime;

%% Main Loop
T = Inputs.Time.StartTime - Inputs.Time.dT;

while T < Inputs.Time.EndTime
    
    % Move to next timestep
    T = T + Inputs.Time.dT;
    Cell.WetLastTimestep = Cell.Wet;
    
    % Update bed composition
    if Inputs.Sed.SedType == 2
        % NEEDS MORE THOUGHT HERE ABOUT MIXING MODEL...
        % Fi = (Fi * DA + Delta_i_tot * DT) / DA; % Could mix new sed in before moving from active to sub-surface??? Timestep sensitive?
        Cell.SubSurfFlux_i = (max(-Cell.Delta_tot,0) * ones(1,Frac.NFracs)) .* Cell.BulkFi - ...
                        (max(Cell.Delta_tot,0) * ones(1,Frac.NFracs)) .* Cell.Fi; % Fractional volumetric flux rate into active layer from subsurface [m3/s/m]
        % Fi = (Fi * DA - SubSurfFlux_i) / DA;
        Cell.Fi = (Cell.Fi * Inputs.Sed.DA + (Cell.Delta_i_tot - Cell.SubSurfFlux_i) * Inputs.Time.dT) / Inputs.Sed.DA;
    end

    % Update bed level
    Cell.Z = Cell.Z + Cell.Delta_tot * Inputs.Time.dT / Inputs.Sed.Porosity;
    
    % Calculate basic hydraulics
    Cell.H = zeros(Cell.NCells,1);
    Cell.U = zeros(Cell.NCells,1);
    Cell.Tau_S = zeros(Cell.NCells,1);
    
    WL = SetQ(Cell.Width,Cell.Z,Cell.WetLastTimestep,Inputs.Hyd);
    Cell.Wet = (WL-Cell.Z >= Inputs.Hyd.DryFlc) | (Cell.WetLastTimestep & (WL-Cell.Z >= Inputs.Hyd.DryFlc / 2));
    Cell.H(Cell.Wet) = WL - Cell.Z(Cell.Wet);
    Cell.U(Cell.Wet) = Cell.H(Cell.Wet).^(2/3) * Inputs.Hyd.Slope^0.5 / Inputs.Hyd.Roughness;
    Cell.Tau_S(Cell.Wet) = Inputs.Hyd.Rho_W * Inputs.Hyd.g * Cell.U(Cell.Wet).^2 * Inputs.Hyd.Roughness^2 ./ Cell.H(Cell.Wet).^(1/3);

    % Calculate secondary flow
    Cell.AlphaSpiral = NaN(Cell.NCells,1);
    Cell.SpiralIntensity = zeros(Cell.NCells,1);
    Cell.Tau_N = zeros(Cell.NCells,1);
    
    Cell.AlphaSpiral(Cell.Wet) = min(sqrt(Inputs.Hyd.g) ./ (Inputs.Hyd.Kappa * (Cell.H(Cell.Wet).^(1/6) / Inputs.Hyd.Roughness)), 0.5); % Equation 9.156 in Delft3D-FLOW User Manual or Eq8 Kalkwijk & Booji 1986
    Cell.SpiralIntensity(Cell.Wet) = Cell.H(Cell.Wet) .* Cell.U(Cell.Wet) / Inputs.Hyd.Radius; % Eq 9.160 (I_be (intensity due to bend) only as I_ce (coriolis) is negligable)
    
    Cell.Tau_N(Cell.Wet) = -2 * Inputs.Opt.Spiral.ESpiral * Inputs.Hyd.Rho_W * Cell.AlphaSpiral(Cell.Wet).^2 .* (1 - Cell.AlphaSpiral(Cell.Wet) / 2) .* Cell.U(Cell.Wet) .* Cell.SpiralIntensity(Cell.Wet); % Eq9.145 D3D or Eq26 Kalkwijk & Booji 1986
    Cell.Tau_Tot = sqrt(Cell.Tau_S.^2 + Cell.Tau_N.^2);
    Cell.UStar = sqrt(Cell.Tau_Tot / Inputs.Hyd.Rho_W);

    % Sediment properties
    Cell.Dg_phi = sum((ones(Cell.NCells,1)*Frac.Di_phi) .* Cell.Fi, 2);
    Cell.Dg_m = 2.^-Cell.Dg_phi ./ 1000;
    Cell.SigmaG_phi = sqrt(sum(Cell.Fi .* ((ones(Cell.NCells,1)*Frac.Di_phi) - (Cell.Dg_phi*ones(1,Frac.NFracs))).^2, 2));

    % Calculate voumetric sediment transport due to flow (in cell centres)
    Cell.Active = (WL - Cell.Z) >= Inputs.Sed.SedThr;
    Edge.Active = [0; Cell.Active(1:end-1) .* Cell.Active(2:end); 0];
    if Inputs.Opt.ST.Formula == 1
        Cell.qsiTot_flow = WilcockCrowe(Inputs.Sed.Rho_S, Inputs.Hyd.Rho_W, Inputs.Hyd.g, Cell.NCells, Frac.NFracs, Frac.Di_m, Frac.SandFrac, Cell.Fi, Cell.Dg_m, Cell.Tau_Tot, Cell.UStar);
    end
    Cell.qsiTot_flow(~Cell.Active) = 0;
    Cell.qsS_flow_kg = zeros(Cell.NCells,1);
    Cell.qsS_flow_kg(Cell.Active) = Inputs.Sed.Rho_S * sum(Cell.qsiTot_flow(Cell.Active), 2) .* (Cell.Tau_S(Cell.Active) ./ Cell.Tau_Tot(Cell.Active)); % Streamwise mass transport rate [kg/s/m]
    Cell.qsiN_flow = Cell.qsiTot_flow .* ((Cell.Tau_N ./ Cell.Tau_Tot) * ones(1, Frac.NFracs));
    Cell.qsiN_flow(isnan(Cell.qsiN_flow)) = 0;

    % Calculate cell edge parameters
    if Inputs.Opt.UpwindBedload % Upwind bedload
        C2E_Weights = (Cell.Tau_N(1:end-1) + Cell.Tau_N(2:end)) > 0;
    else % Central
        C2E_Weights = 0.5 * ones(Cell.NCells - 1, 1);
    end
    
    Edge.qsiTot_flow = Centre2Edge(Cell.qsiTot_flow, C2E_Weights);
    Edge.qsiTot_flow(~Edge.Active,:) = 0;
    Edge.qsiN_flow = Centre2Edge(Cell.qsiN_flow, C2E_Weights);
    Edge.qsiN_flow(~Edge.Active,:) = 0;
    Edge.H = Centre2Edge(Cell.H, C2E_Weights);
    Edge.Tau_Tot = Centre2Edge(Cell.Tau_Tot, C2E_Weights);
    Edge.Dg_m = Centre2Edge(Cell.Dg_m, C2E_Weights);
    
    % Calculate sediment transport due to bed slope
    Edge.Slope = [0; (Cell.Z(2:end) - Cell.Z(1:end-1)) ./ (Cell.N(2:end) - Cell.N(1:end-1)); 0];
    % Note: at the moment, for the purposes of slope analysis we are
    %       assuming unadjusted bedload transport vector is normal to
    %       cross-section and that downstream slope << transverse slope
    %       (i.e. secondary flow neglected, downstream slope neglected for 
    %       bed slope calc). 
    switch Inputs.Opt.Slope.Formula
        case 0
            % No bed slope influence on transport
            Edge.qsiN_slope = zeros(Edge.NEdges, Frac.NFracs);
        case 1
            % Talmon et al 1995
            Edge.ShieldsStressi = (Edge.Tau_Tot * ones(1, Frac.NFracs)) ./ ((Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * Inputs.Hyd.g * (ones(Edge.NEdges, 1) * Frac.Di_m));
            Edge.fTheta = Inputs.Opt.Slope.A_sh * Edge.ShieldsStressi.^Inputs.Opt.Slope.B_sh .* ...
                     (((ones(Edge.NEdges, 1) * Frac.Di_m)) ./ (Edge.H * ones(1, Frac.NFracs))).^Inputs.Opt.Slope.C_sh .* ...
                     ((ones(Edge.NEdges,1) * Frac.Di_m) ./ (Edge.Dg_m * ones(1,Frac.NFracs))).^Inputs.Opt.Slope.D_sh;
            Edge.qsiN_slope = - Edge.qsiTot_flow .* (1 ./ Edge.fTheta) .* (Edge.Slope * ones(1, Frac.NFracs));
            Edge.qsiN_slope(isnan(Edge.qsiN_slope)) = 0;
    end
    
    % Erosion/deposition due to flow and associated effects (i.e. everything except bank erosion)
    Cell.Delta_i_flow = Edge.qsiN_flow(1:end-1,:) - Edge.qsiN_flow(2:end,:); % Fractional volumetric flux rate into cell from neighboring cells due to flow [m3/s/m]
    Cell.Delta_i_slope = Edge.qsiN_slope(1:end-1,:) - Edge.qsiN_slope(2:end,:);
    
    % Calculate sediment flux rate due to bank erosion
    % Identify banks
    Edge.IsBank = IdentifyBanks(Inputs.Opt.Bank.ID, Cell, Edge);
    
    % Apply bank erosion stencil
    Bank = BankStencil(Inputs.Opt.Bank.Stencil, Cell, Edge);
    
    % Calculate flux
    Cell.Delta_i_bank = BankFlux(Inputs.Opt.Bank.Flux, Cell, Edge, Frac, Inputs.Time.dT, Bank);
    Cell.Delta_i_tot = Cell.Delta_i_flow + Cell.Delta_i_slope + Cell.Delta_i_bank;
    Cell.Delta_tot = sum(Cell.Delta_i_tot, 2); % total volumetric flux rate into cell from neighboring cells [m3/s/m]
    
    % Outputs
    if T >= PlotT + Inputs.Outputs.PlotInt
        PlotT = PlotT + Inputs.Outputs.PlotInt;
        UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T)
    end
    if T >= DiagT + Inputs.Outputs.DiagInt
        DiagT = DiagT + Inputs.Outputs.DiagInt;
        fprintf('T=%gs, Q=%gm^3/s, q_s=%.2em3/s\n', T, Inputs.Hyd.Flow, sum(Cell.qsS_flow_kg))
    end
end

%% Final outputs
UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T)
