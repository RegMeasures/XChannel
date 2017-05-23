function [Cell, Edge, Frac, Bank] = InitialiseVariables(Inputs)
%INITIALISEVARIABLES   Create variables ready for XChannel loops to run
%Create/populate structs and fields where necessary before XChannel model
%starts iterating. Not all variables are created - only ones which are
%needed prior to their definition elsewhere - e.g. variables used for plot 
%generation or calculation inputs.
%
%   [Cell, Edge, Frac, Bank] = INITIALISEVARIABLES(Inputs)
%
%   Inputs:
%      Inputs = Structure array created by ReadModelInputs
%
%   Outputs:
%      Cell = Structure array to hold variables of cell center properties.
%             The initial fields created to by INITIALISEVARIABLES are:
%         .NCells   = number of cells across cross-section
%         .N        = cell center location, cross-channel dir [m]
%                     (NCells x 1)
%         .Width    = cell width [m] (NCells x 1)
%         .Z        = cell bed elevation [m] (NCells x 1)
%         .Zinitial = cell initial bed elevation [m] (NCells x 1)
%         .Zfinal   = observed final bed elevation if supplied [m]
%                     (NCells x 1)
%         .EroStore = Stored bank erosion (only if StoredBE = 1) [m]
%                     (NCells x 1)
%         .Wet      = boolean flag for wet cells (NCells x 1)
%         .U        = depth averaged streamwise velocity in cell [m/s]
%                     (NCells x 1) 
%         .Tau_S    = streamwise bed shear stress in cell [N/m2] 
%                     (NCells x 1)
%         .Fi       = proportion of each sediment fraction in acive layer
%                     of cell (NCells x NFracs)
%         .BulkFi   = proportion of each sediment fraction in subsurface
%                     of cell (NCells x NFracs)
%         .SubDg_m  = geometric mean grain size of subsurface [m](NCellsx1) 
%         .Dg_phi   = geometric mean grain size of active layer [phi]
%                     (NCells x 1) 
%         .Dg_m     = geometric mean grain size of active layer [m]
%                     (NCells x 1) 
%         .qsS_flow_kg = Streamwise bedload transport rate [kg/s/m] 
%                     (NCells x 1) 
%         .Delta_tot = Rate of deposition in each cell [m3/m/s](NCells x 1) 
%         .Delta_i_tot = Fractional rate of deposition in each cell
%                     [m3/m/s] (NCells x NFracs) 
%      Edge = Structure array to hold variables of cell edge properties
%         .NEdges   = number of edges across cross-section = NCells + 1
%         .N        = cell edge location, cross-channel dir (NEdges x 1)
%         .IsBank   = flag for whether edge is a bank (NEdges x 1)
%      Frac = Structure array of sediment fraction properties
%         .NFracs   = number of sediment fractions
%         .Di_m     = size of sediment fraction [m] (1 x NFracs)
%         .Di_phi   = size of sediment fraction [phi] (1 x NFracs)
%         .SandFrac = flag for fractions which are 'sand' for purposes of
%                     Wilcock-Crowe bedload transport (1 x NFracs)
%      Bank = Structure array of bank properties
%         .NBanks   = number of identified banks
%         .Top      = bank top cell number (NBanks x 1)
%         .Bottom   = bank bottom cell number (NBanks x 1)
%
%   Notes: 
%      Several fields are initialised with zero or dummy values.
%      Additional fields not initialised here are added to Cell, Edge, and
%      Bank structure arrays by other parts of XChannel.
%
%   See also: XCHANNEL, READMODELINPUTS.

%% Initialise geometry
Cell.NCells = size(Inputs.Hyd.InitialGeometry,1);
Edge.NEdges = Cell.NCells + 1;
Cell.N = Inputs.Hyd.InitialGeometry(:,1);
Edge.N = [Cell.N(1,1) - (Cell.N(2,1) - Cell.N(1,1)) / 2;
               (Cell.N(1:end-1,1) + Cell.N(2:end,1)) / 2;
               Cell.N(end,1) + (Cell.N(end,1) - Cell.N(end-1,1)) / 2];
Cell.Width = Edge.N(2:end,1) - Edge.N(1:end-1,1);
Cell.Z = Inputs.Hyd.InitialGeometry(:,2);
Cell.Zinitial = Inputs.Hyd.InitialGeometry(:,2);
if size(Inputs.Hyd.InitialGeometry,2) >= 3
    Cell.Zfinal = Inputs.Hyd.InitialGeometry(:,3);
else
    Cell.Zfinal = NaN(Cells.NCells,1);
end
if Inputs.Bank.Update.StoredBE % If storing bank erosion (Nicholas 2013 approach)
    Cell.EroStore = zeros(Cell.NCells,1);
end

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
    if sum(Inputs.Sed.SedSize(:,2))>1.1 || sum(Inputs.Sed.SedSize(:,2))<0.9
        error('Surface grain size distribution must sum to 1')
    end
    Cell.Fi = ones(Cell.NCells,1) * Inputs.Sed.SedSize(:,2)' / sum(Inputs.Sed.SedSize(:,2));
    if sum(Inputs.Sed.SedSize(:,3))>1.1 || sum(Inputs.Sed.SedSize(:,3))<0.9
        error('Sub-surface grain size distribution must sum to 1')
    end
    Cell.BulkFi = ones(Cell.NCells,1) * Inputs.Sed.SedSize(:,3)' / sum(Inputs.Sed.SedSize(:,3));
else
    error('SedSize input type not valid')
end
Frac.NFracs = size(Frac.Di_m,2);
Frac.Di_phi = -log2(Frac.Di_m * 1000);
if Inputs.ST.Formula == 2 % Which fractions are sand - required for Wilcock-Crowe bedload formula
    Frac.SandFrac = Frac.Di_m <= 0.002;
end
Cell.SubDg_m = 2.^-sum((ones(Cell.NCells,1)*Frac.Di_phi) .* Cell.BulkFi, 2) ./1000;
Cell.Dg_phi = sum((ones(Cell.NCells,1)*Frac.Di_phi) .* Cell.Fi, 2);
Cell.Dg_m = 2.^-Cell.Dg_phi ./ 1000;

%% Initialise transport
Cell.qsS_flow_kg = NaN(Cell.NCells,1);

%% Initialise banks
Edge.IsBank = zeros(Edge.NEdges,1);
Bank.Top = [];
Bank.Bottom = [];
Bank.NBanks = 0;

%% Initialise morpho
Cell.Delta_tot = zeros(Cell.NCells,1); % total volumetric flux rate into cell from neighboring cells [m3/s/m]
Cell.Delta_i_tot = zeros(Cell.NCells,Frac.NFracs); % fractional volumetric flux rate into cell from neighboring cells [m3/s/m]

%% Initialise Time


end

