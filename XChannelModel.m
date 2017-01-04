function [FinalXS, WL] = XChannelModel(Inputs)
%XCHANNELMODEL Single cross-section 1D (cross-channel) morphological model
%Morphological model for simulating the evolution of cross-section shape
%for a single cross section
%
%   [FinalXS] = XChannelModel() Runs model prompting user to select 
%   input file
%
%   [FinalXS] = XChannelModel(FileName) Runs model with inputs from
%   Filename, where filename is a text file in standard format
%
%   [FinalXS] = XChannelModel(Inputs) Runs model with inputs contained in
%   struct created by ReadModelInputs
%
%   FinalXS is the final cross-section profile at the end of the simulation.
%
%   See also: READMODELINPUTS

%% Program setup
if ~isdeployed
    addpath('Functions')
end

%% Read model data and options from model input file if required
if ~exist('Inputs','var')
        [Inputs] = ReadModelInputs;
else
    if ~isstruct(Inputs)
        [Inputs] = ReadModelInputs(Inputs);
    end
end

%% Initialise the main variable structs
[Cell, Edge, Frac, Bank] = InitialiseVariables(Inputs);

%% Set up cross-section plot
PlotSed = (Inputs.Sed.SedType == 2);
if Inputs.Outputs.PlotInt > 0
    XsFigure = PlotXS(Cell,Edge,Bank,NaN,0,Inputs.Hyd.Flow,PlotSed);
end

%% Open file for video output
if Inputs.Outputs.VideoOut == 1 && Inputs.Outputs.PlotInt > 0;
    vidObj = VideoWriter(Inputs.FileName(1:end-4),'MPEG-4');
    open(vidObj);
end

%% Create folder for snapshots
if Inputs.Outputs.CsvInt > 0;
    [SnapshotDir,ScenarioTitle,~] = fileparts(Inputs.FileName);
    SnapshotDir = [SnapshotDir,'\snapshots_',ScenarioTitle];
    if ~exist(SnapshotDir,'dir')
        mkdir(SnapshotDir)
    end
end

%% Initialise Time
PlotT = Inputs.Time.StartTime;
DiagT = Inputs.Time.StartTime;
CsvT  = Inputs.Time.StartTime;

T = Inputs.Time.StartTime - Inputs.Time.dT;

%% Main Loop
while T < Inputs.Time.EndTime  
    %% Move to next timestep
    T = T + Inputs.Time.dT;
    Cell.WetLastTimestep = Cell.Wet;
    
    %% Update bed level
    Cell.Z = Cell.Z + Cell.Delta_tot * Inputs.Time.dT ./ Cell.Width / (1-Inputs.Sed.Porosity);  
    
    %% Update bed composition
    if Inputs.Sed.SedType == 2 % if graded sediment
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
        Cell.SigmaG_phi = sqrt(sum(Cell.Fi .* ...
                                   ((ones(Cell.NCells,1)*Frac.Di_phi) - ...
                                    (Cell.Dg_phi*ones(1,Frac.NFracs))).^2, 2));
    end
    
    %% Get new flow (if time series)
    if Inputs.Hyd.FlowType == 2;
        Inputs.Hyd.Flow = interp1(Inputs.Hyd.FlowTS(:,1), ...
                                  Inputs.Hyd.FlowTS(:,2),T);
    end
    
    %% Calculate basic hydraulics

    WL = SetQ(Cell.Width,Cell.Z,Cell.WetLastTimestep,Inputs.Hyd);
    
    [Cell.Chezy, Cell.H, Cell.U, Cell.Tau_S, Cell.Wet] = ...
        BasicHydraulics(Cell.Z, Cell.WetLastTimestep, Inputs.Hyd, WL);
    
    %% Calculate transverse shear stress due to secondary flow
    [Cell.Tau_N] =  SecondaryFlow(Inputs.Hyd, Cell);
    
    Cell.Tau_Tot = sqrt(Cell.Tau_S.^2 + Cell.Tau_N.^2);
    Cell.UStar = sqrt(Cell.Tau_Tot / Inputs.Hyd.Rho_W);
    
    %% Calculate fractional volumetric bedload transport due to flow
    Cell.ShieldsStress_i = (Cell.Tau_Tot * ones(1,Frac.NFracs)) ./ ...
                           ((Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * ...
                            Inputs.Hyd.g * ...
                            (ones(Cell.NCells,1) * Frac.Di_m));
                     
    Cell.Active = (WL - Cell.Z) >= Inputs.Sed.SedThr; 
    % Cells are only active for transport if depth is greater than SedThr
    
    Edge.Active = [0; Cell.Active(1:end-1) .* Cell.Active(2:end); 0]; 
    % Edge is active for transport only if cells on both sides are active
    
    [Cell.qsiTot_flow,Cell.ThetaCrit_i] = BedLoad(Inputs, Cell, Frac);
    
    Cell.qsTot_flow = sum(Cell.qsiTot_flow, 2);
    
    % qsS_flow_kg = Streamwise mass transport rate [kg/s/m]
    Cell.qsS_flow_kg = zeros(Cell.NCells,1);
    Cell.qsS_flow_kg(Cell.Active) = ...
        Inputs.Sed.Rho_S *  Cell.qsTot_flow(Cell.Active,:) .* ...
        (Cell.Tau_S(Cell.Active) ./ Cell.Tau_Tot(Cell.Active)); 
    
    Cell.qsiN_flow = Cell.qsiTot_flow .* ...
                     ((Cell.Tau_N ./ Cell.Tau_Tot) * ones(1, Frac.NFracs));
    Cell.qsiN_flow(isnan(Cell.qsiN_flow)) = 0;

    %% Calculate cell edge parameters
    % note: the level of upwinding vs central is controlled by 
    % Inputs.ST.UpwindBedload
    
    % Identify upwind cells
    C2E_Weights = (Cell.Tau_N(1:end-1) + Cell.Tau_N(2:end)) > 0;
    % Set weightings for LH cells (RH weightings = 1-C2E_Weights)
    C2E_Weights = Inputs.ST.UpwindBedload * C2E_Weights + ...
                  (1-Inputs.ST.UpwindBedload) * (C2E_Weights==0);
    
    % Calculate cell edge parameters using pre-calculated weightings
    Edge.H               = Centre2Edge(Cell.H, C2E_Weights, Edge.Active);
    Edge.Tau_Tot         = Centre2Edge(Cell.Tau_Tot, C2E_Weights, ...
                                       Edge.Active);    
    Edge.ShieldsStress_i = Centre2Edge(Cell.ShieldsStress_i, ...
                                       C2E_Weights, Edge.Active);
    Edge.qsiTot_flow     = Centre2Edge(Cell.qsiTot_flow, C2E_Weights, ...
                                       Edge.Active);
    Edge.qsiN_flow       = Centre2Edge(Cell.qsiN_flow, C2E_Weights, ...
                                       Edge.Active);
    Edge.ThetaCrit_i     = Centre2Edge(Cell.ThetaCrit_i, C2E_Weights, ...
                                       Edge.Active);
    Edge.Dg_m            = Centre2Edge(Cell.Dg_m, C2E_Weights, ...
                                       Edge.Active);
    
    %% Calculate transverse sediment transport due to bed slope
    Edge.Slope = [0; (Cell.Z(2:end) - Cell.Z(1:end-1)) ./ ...
                  (Cell.N(2:end) - Cell.N(1:end-1)); 0];
    Edge.qsiN_slope = BedSlope(Inputs.Slope, Edge, Frac);
    
    %% Fractional erosion/deposition due to flow and associated effects (i.e. everything except bank erosion)
    Cell.Delta_i_flow = Edge.qsiN_flow(1:end-1,:) - Edge.qsiN_flow(2:end,:); % Fractional volumetric flux rate into cell from neighboring cells due to flow [m3/s/m]
    Cell.Delta_i_slope = Edge.qsiN_slope(1:end-1,:) - Edge.qsiN_slope(2:end,:);
    Cell.Delta_flow = sum(Cell.Delta_i_flow, 2); % Total (i.e. all fractions) flux into cell from neighbouring cells due to flow
    Cell.Delta_slope = sum(Cell.Delta_i_slope, 2);
    
    %% Calculate sediment flux rate due to bank erosion
    % Identify banks
    Edge.IsBank = IdentifyBanks(Inputs.Bank.ID, Cell, Edge);
    
    % Apply bank erosion stencil
    Bank = BankStencil(Inputs.Bank.Stencil, Cell, Edge);
    
    % Apply active bank trigger
    Bank.Active = TriggerBanks(Inputs.Bank.Trigger, Cell, Bank);
    
    % Calculate bank erosion flux
    Cell.Delta_i_bank = BankFlux(Inputs.Bank.Flux, Cell, Frac, ...
                                 Inputs.Time.dT, Bank);
    Cell.Delta_bank = sum(Cell.Delta_i_bank, 2);
    
    %% Calculate total erosion/deposition
    Cell.Delta_i_tot = Cell.Delta_i_flow + Cell.Delta_i_slope + Cell.Delta_i_bank;
    if Inputs.Bank.Update.StoredBE
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
        UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T, Inputs.Hyd.Flow, PlotSed)
        % Write video frame
        if Inputs.Outputs.VideoOut == 1
            if PlotSed
                Frame = getframe(XsFigure.FigureH,[55, 393, 906, 382]);
            else
                Frame = getframe(XsFigure.FigureH,[55, 193, 906, 382]);
            end
            writeVideo(vidObj,Frame);
        end
    end
    
    % Diagnistics output
    if T >= DiagT + Inputs.Outputs.DiagInt && Inputs.Outputs.DiagInt > 0
        DiagT = DiagT + Inputs.Outputs.DiagInt;
        fprintf('T=%gs, Q=%gkg/s, q_s=%.2em3/s\n', T, Inputs.Hyd.Flow, sum(Cell.qsS_flow_kg))
    end
    
    % CSV output
    if T >= CsvT + Inputs.Outputs.CsvInt && Inputs.Outputs.CsvInt > 0
        CsvT = CsvT + Inputs.Outputs.CsvInt;
        csvwrite(sprintf('%s\\BedSnapshot_T=%i.out',SnapshotDir,T),[Cell.N,Cell.Z]);
    end
end

%% Final tidying up
% UpdateXsPlot(XsFigure, Cell, Edge, Bank, WL, T, Inputs.Hyd.Flow)

% Save and close figure
if Inputs.Outputs.PlotInt > 0
    saveas(XsFigure.FigureH,Inputs.FileName(1:end-4),'png')
    close(XsFigure.FigureH);
    % Close video file
    if Inputs.Outputs.VideoOut
        close(vidObj)
    end
end

%% Return final bed level to allow fit analysis
FinalXS = [Cell.N,Cell.Z];

if Inputs.Outputs.CsvInt > 0
    csvwrite(sprintf('%s\\BedSnapshot_T=%i.out',SnapshotDir,T),FinalXS);
end

end

