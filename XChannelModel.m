function [FinalXS] = XChannelModel(Inputs)
% Single cross-section 1D (cross-channel) morphological model
%
% [FinalXS] = XChannelModel
%     Run model prompting user to select input file
%
% [FinalXS] = XChannelModel(FileName)
%     Run model with inputs from Filename where filename is a text file in
%     standard format
%
% [FinalXS] = XChannelModel(Inputs)
%     Run model with inputs contained in struct created by ReadModelInputs
%
% FinalXS is the final cross-section profile at the end of the simulation.
%
% FileName is optional name(and path) of input file
% If FileName not specified user will be prompted to select file
% Richard Measures 2015

%% Program setup
if ~isdeployed
    addpath('Functions')
end

%% Read input file if required
if ~isstruct(Inputs)
    FileName = Inputs;
    % Get file name and path if not specified as input
    if ~exist('FileName','var')
        [FileName,FilePath] = uigetfile('*.txt','Select the model input file');
        if isequal(FileName,0)
            error('User selected Cancel')
        end
        FileName = fullfile(FilePath,FileName);
    end

    % Read model data and options from intput file
    [Inputs] = ReadModelInputs(FileName);
end

%% Initialise the main variable structs
[Cell, Edge, Frac, Bank] = InitialiseVariables(Inputs);

%% Set up cross-section plot
if Inputs.Outputs.PlotInt > 0
    XsFigure = PlotXS(Cell,Edge,Bank,NaN,0,Inputs.Hyd.Flow);
end    
PlotT = Inputs.Time.StartTime;
DiagT = Inputs.Time.StartTime;


%% Open file for video output
if Inputs.Outputs.VideoOut == 1 && Inputs.Outputs.PlotInt > 0;
    vidObj = VideoWriter(Inputs.FileName(1:end-4),'MPEG-4');
    open(vidObj);
end

%% Main Loop
T = Inputs.Time.StartTime - Inputs.Time.dT;

while T < Inputs.Time.EndTime  
    %% Move to next timestep
    T = T + Inputs.Time.dT;
    Cell.WetLastTimestep = Cell.Wet;
    
    %% Update bed level
    Cell.Z = Cell.Z + Cell.Delta_tot * Inputs.Time.dT ./ Cell.Width / Inputs.Sed.Porosity;  
    
    %% Update bed composition
    if Inputs.Sed.SedType == 2
        % NEEDS MORE THOUGHT HERE ABOUT MIXING MODEL...
        % Fi = (Fi * DA + Delta_i_tot * DT) / DA; % Could mix new sed in before moving from active to sub-surface??? Timestep sensitive?
        Cell.SubSurfFlux_i = (max(-Cell.Delta_tot,0) * ones(1,Frac.NFracs)) .* Cell.BulkFi - ...
                        (max(Cell.Delta_tot,0) * ones(1,Frac.NFracs)) .* Cell.Fi; % Fractional volumetric flux rate into active layer from subsurface [m3/s/m]
        % Fi = (Fi * DA - SubSurfFlux_i) / DA;
        Cell.Fi = (Cell.Fi .* (Cell.Width * ones(1,Frac.NFracs)) * Inputs.Sed.DA + ...
                   (Cell.Delta_i_tot + Cell.SubSurfFlux_i) * Inputs.Time.dT ) ./ ...
                  (Inputs.Sed.DA * (Cell.Width * ones(1,Frac.NFracs)));
        Cell.Fi(Cell.Fi<0) = 0;
        Cell.Fi = Cell.Fi ./ (sum(Cell.Fi,2) * ones(1,Frac.NFracs)); % To prevent small errors propogating?
    end

    %% Update sediment properties (if graded sediment)
    if Inputs.Sed.SedType == 2
        Cell.Dg_phi = sum((ones(Cell.NCells,1)*Frac.Di_phi) .* Cell.Fi, 2);
        Cell.Dg_m = 2.^-Cell.Dg_phi ./ 1000;
        Cell.SigmaG_phi = sqrt(sum(Cell.Fi .* ((ones(Cell.NCells,1)*Frac.Di_phi) - (Cell.Dg_phi*ones(1,Frac.NFracs))).^2, 2));
    end
    
    %% Get new flow (if time series)
    if Inputs.Hyd.FlowType == 2;
        Inputs.Hyd.Flow = interp1(Inputs.Hyd.FlowTS(:,1),Inputs.Hyd.FlowTS(:,2),T);
    end
    
    %% Calculate basic hydraulics
    Cell.H = zeros(Cell.NCells,1);
    Cell.U = zeros(Cell.NCells,1);
    Cell.Tau_S = zeros(Cell.NCells,1);
    
    WL = SetQ(Cell.Width,Cell.Z,Cell.WetLastTimestep,Inputs.Hyd);
    Cell.Wet = (WL-Cell.Z >= Inputs.Hyd.DryFlc) | (Cell.WetLastTimestep & (WL-Cell.Z >= Inputs.Hyd.DryFlc / 2));
    Cell.H(Cell.Wet) = WL - Cell.Z(Cell.Wet);
    Cell.U(Cell.Wet) = Cell.H(Cell.Wet).^(2/3) * Inputs.Hyd.Slope^0.5 / Inputs.Hyd.Roughness;
    Cell.Tau_S(Cell.Wet) = Inputs.Hyd.Rho_W * Inputs.Hyd.g * Cell.U(Cell.Wet).^2 * Inputs.Hyd.Roughness^2 ./ Cell.H(Cell.Wet).^(1/3);

    %% Calculate secondary flow
    Cell.AlphaSpiral = NaN(Cell.NCells,1);
    Cell.SpiralIntensity = zeros(Cell.NCells,1);
    Cell.Tau_N = zeros(Cell.NCells,1);
    
    Cell.AlphaSpiral(Cell.Wet) = min(sqrt(Inputs.Hyd.g) ./ (Inputs.Hyd.Kappa * (Cell.H(Cell.Wet).^(1/6) / Inputs.Hyd.Roughness)), 0.5); % Equation 9.156 in Delft3D-FLOW User Manual or Eq8 Kalkwijk & Booji 1986
    Cell.SpiralIntensity(Cell.Wet) = Cell.H(Cell.Wet) .* Cell.U(Cell.Wet) / Inputs.Hyd.Radius; % Eq 9.160 (I_be (intensity due to bend) only as I_ce (coriolis) is negligable)
    
    Cell.Tau_N(Cell.Wet) = -2 * Inputs.Opt.Spiral.ESpiral * Inputs.Hyd.Rho_W * Cell.AlphaSpiral(Cell.Wet).^2 .* (1 - Cell.AlphaSpiral(Cell.Wet) / 2) .* Cell.U(Cell.Wet) .* Cell.SpiralIntensity(Cell.Wet); % Eq9.145 D3D or Eq26 Kalkwijk & Booji 1986
    Cell.Tau_Tot = sqrt(Cell.Tau_S.^2 + Cell.Tau_N.^2);
    Cell.UStar = sqrt(Cell.Tau_Tot / Inputs.Hyd.Rho_W);
    
    %% Calculate fractional volumetric bedload transport due to flow (in cell centres)
    Cell.ShieldsStress_i = (Cell.Tau_Tot * ones(1,Frac.NFracs)) ./ ...
                          ((Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * Inputs.Hyd.g * (ones(Cell.NCells,1)*Frac.Di_m));
                     
    Cell.Active = (WL - Cell.Z) >= Inputs.Sed.SedThr;
    Edge.Active = [0; Cell.Active(1:end-1) .* Cell.Active(2:end); 0];
    
    [Cell.qsiTot_flow,Cell.ThetaCrit_i] = BedLoad(Inputs, Cell, Frac);
    
    Cell.qsTot_flow = sum(Cell.qsiTot_flow, 2);
    Cell.qsS_flow_kg = zeros(Cell.NCells,1);
    Cell.qsS_flow_kg(Cell.Active) = Inputs.Sed.Rho_S *  Cell.qsTot_flow(Cell.Active,:) .* (Cell.Tau_S(Cell.Active) ./ Cell.Tau_Tot(Cell.Active)); % Streamwise mass transport rate [kg/s/m]
    Cell.qsiN_flow = Cell.qsiTot_flow .* ((Cell.Tau_N ./ Cell.Tau_Tot) * ones(1, Frac.NFracs));
    Cell.qsiN_flow(isnan(Cell.qsiN_flow)) = 0;

    %% Calculate cell edge parameters
    if Inputs.Opt.UpwindBedload % Upwind bedload
        C2E_Weights = (Cell.Tau_N(1:end-1) + Cell.Tau_N(2:end)) > 0;
    else % Central
        C2E_Weights = 0.5 * ones(Cell.NCells - 1, 1);
    end
    
    Edge.H = Centre2Edge(Cell.H, C2E_Weights);
    Edge.Tau_Tot = Centre2Edge(Cell.Tau_Tot, C2E_Weights);    
    Edge.ShieldsStress_i = Centre2Edge(Cell.ShieldsStress_i, C2E_Weights);
    Edge.qsiTot_flow = Centre2Edge(Cell.qsiTot_flow, C2E_Weights);
    Edge.qsiTot_flow(~Edge.Active,:) = 0;
    Edge.qsiN_flow = Centre2Edge(Cell.qsiN_flow, C2E_Weights);
    Edge.qsiN_flow(~Edge.Active,:) = 0;
    Edge.ThetaCrit_i = Centre2Edge(Cell.ThetaCrit_i, C2E_Weights);
    Edge.Dg_m = Centre2Edge(Cell.Dg_m, C2E_Weights);
    
    %% Calculate sediment transport due to bed slope
    Edge.Slope = [0; (Cell.Z(2:end) - Cell.Z(1:end-1)) ./ (Cell.N(2:end) - Cell.N(1:end-1)); 0];
    Edge.qsiN_slope = BedSlope(Inputs, Edge, Frac);
    
    %% Fractional erosion/deposition due to flow and associated effects (i.e. everything except bank erosion)
    Cell.Delta_i_flow = Edge.qsiN_flow(1:end-1,:) - Edge.qsiN_flow(2:end,:); % Fractional volumetric flux rate into cell from neighboring cells due to flow [m3/s/m]
    Cell.Delta_i_slope = Edge.qsiN_slope(1:end-1,:) - Edge.qsiN_slope(2:end,:);
    Cell.Delta_flow = sum(Cell.Delta_i_flow, 2); % Total (i.e. all fractions) flux into cell from neighbouring cells due to flow
    Cell.Delta_slope = sum(Cell.Delta_i_slope, 2);
    
    %% Calculate sediment flux rate due to bank erosion
    % Identify banks
    Edge.IsBank = IdentifyBanks(Inputs.Opt.Bank.ID, Cell, Edge);
    
    % Apply bank erosion stencil
    Bank = BankStencil(Inputs.Opt.Bank.Stencil, Cell, Edge);
    
    % Calculate bank erosion flux
    Cell.Delta_i_bank = BankFlux(Inputs.Opt.Bank.Flux, Cell, Edge, Frac, Inputs.Time.dT, Bank);
    Cell.Delta_bank = sum(Cell.Delta_i_bank, 2);
    
    %% Calculate total erosion/deposition
    Cell.Delta_i_tot = Cell.Delta_i_flow + Cell.Delta_i_slope + Cell.Delta_i_bank;
    if Inputs.Opt.Bank.Flux.StoredBE
        % Store bank erosion if intermittent update option selected
        Cell.Delta_store = StoreErosion(Inputs, Cell, Bank);
        Cell.EroStore = Cell.EroStore - Cell.Delta_store * Inputs.Time.dT;
        Cell.Delta_tot = Cell.Delta_flow + Cell.Delta_slope + Cell.Delta_bank + Cell.Delta_store;
    else
        Cell.Delta_tot = Cell.Delta_flow + Cell.Delta_slope + Cell.Delta_bank;
    end
    
    %% Outputs
    
    if T >= PlotT + Inputs.Outputs.PlotInt && Inputs.Outputs.PlotInt > 0
        % Update plot
        PlotT = PlotT + Inputs.Outputs.PlotInt;
        UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T, Inputs.Hyd.Flow)
        % Write video frame
        if Inputs.Outputs.VideoOut == 1
            Frame = getframe(XsFigure.FigureH,[55, 393, 906, 382]);
            writeVideo(vidObj,Frame);
        end
    end
    
    % Diagnistics output
    if T >= DiagT + Inputs.Outputs.DiagInt && Inputs.Outputs.DiagInt > 0
        DiagT = DiagT + Inputs.Outputs.DiagInt;
        fprintf('T=%gs, Q=%gm^3/s, q_s=%.2em3/s\n', T, Inputs.Hyd.Flow, sum(Cell.qsS_flow_kg))
    end
end

%% Final tidying up
% UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T, Inputs.Hyd.Flow)

% Close figure
if Inputs.Outputs.PlotInt > 0
    close(XsFigure.FigureH);
    % Close video file
    if Inputs.Outputs.VideoOut == 1
        close(vidObj)
    end
end

%% Return final bed level to allow fit analysis
FinalXS = [Cell.N,Cell.Z];

end

