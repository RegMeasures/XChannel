function [Cell, Edge, Frac, Bank] = InitialiseVariables(Inputs)
% Create variable structures ready for model to run
% Populate/calculate variables where possible/necessary before model
% starts

% note: Not all variables are created - only ones which are needed for
% plotting or calculation

%% Initialise geometry
Cell.NCells = size(Inputs.Hyd.InitialGeometry,1);
Edge.NEdges = Cell.NCells + 1;
Cell.N = Inputs.Hyd.InitialGeometry(:,1);
Edge.N = [Cell.N(1) - (Cell.N(2) - Cell.N(1)) / 2;
               (Cell.N(1:end-1) + Cell.N(2:end)) / 2;
               Cell.N(end) + (Cell.N(end) - Cell.N(end-1)) / 2];
Cell.Width = Edge.N(2:end) - Edge.N(1:end-1);
Cell.Z = Inputs.Hyd.InitialGeometry(:,2);

%% Initialise hydraulics
Cell.Wet = ones(Cell.NCells,1);
Cell.U = NaN(Cell.NCells,1);
Cell.Tau_S = NaN(Cell.NCells,1);

%% Initialise sediment
if Inputs.Sed.SedType == 1 % uniform sediment
    Frac.Di_m = Inputs.Sed.SedSize;
    Cell.Fi = ones(Cell.NCells,1);
    Cell.BulkFi = Cell.Fi;
elseif Inputs.Sed.SedType == 2 % graded sediment
    Frac.Di_m = Inputs.Sed.SedSize(:,1)';
    Cell.Fi = ones(Cell.NCells,1) * Inputs.Sed.SedSize(:,2)';
    Cell.BulkFi = ones(Cell.NCells,1) * Inputs.Sed.SedSize(:,3)';
end
Frac.NFracs = size(Frac.Di_m,2);
Frac.Di_phi = -log2(Frac.Di_m * 1000);
if Inputs.Opt.ST.Formula == 1
    Frac.SandFrac = Frac.Di_m <= 0.002;
end
Cell.SubDg_m = 2.^-sum((ones(Cell.NCells,1)*Frac.Di_phi) .* Cell.BulkFi, 2) ./1000;
Cell.Dg_m = 2.^-sum((ones(Cell.NCells,1)*Frac.Di_phi) .* Cell.Fi, 2) ./1000;

%% Initialise transport
Cell.qsS_flow_kg = NaN(Cell.NCells,1);

%% Initialise banks
Edge.IsBank = zeros(Edge.NEdges,1);
Bank.Top = [];
Bank.Bottom = [];
Bank.NBanks = 0;

%% Initialise morpho
Cell.Delta_tot = zeros(Cell.NCells,1);
Cell.Delta_i_tot = zeros(Cell.NCells,Frac.NFracs);

%% Initialise Time


end

