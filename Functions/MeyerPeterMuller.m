function [qsiTot_flow, ThetaCrit_i] = MeyerPeterMuller(Inputs, Cell, ...
                                                       Frac, HidExp)
%MEYERPETERMULLER   Bedload transport rate using Meyer-Peter-Muller (1948)
%Fractional volumetric transport rate per unit width calculated using the
%Meyer-Peter-Muller (1948) formulation with a user specified coefficient
%and exponent and externally calculated hiding/exposure correction to
%critical shields stress. Wong and Parker (2006) is a useful reference on
%the application of this formula and the selection of appropriate
%coefficients.
%
%   [qsiTot_flow, ThetaCrit_i] = MEYERPETERMULLER(Inputs, Cell, Frac, ...
%                                                 HidExp)
%   
%   Inputs:
%      Inputs = Struct of model inputs read in by ReadModelInputs.
%               Fields used by MEYERPETERMULLER are:
%         .ST.ThetaCrit    = Critical shields stress
%         .ST.MPMcoef      = MPM coeffieint
%         .ST.MPMexponent  = MPM exponent
%         .Hyd.g           = Acceleration due to gravity [m/s2]
%         .Sed.Rho_S       = Sediment density [kg/m3]
%         .Hyd.Rho_W       = Water density [kg/m3]
%      Cell   = Struct of cell center properties including:
%         .NCells          = number of cells in cross-section
%         .Fi              = NCells x NFrac matrix of proportion of active
%                            layer made up of each fraction, in each cell.
%         .ShieldsStress_i = NCells x NFracs matrix of shields stress
%                            affecting each fraction, in each cell.
%      Frac   = Struct of sediment fraction properties including:
%         .Di_m            = 1 x NFracs matrix of fraction sediment size[m]
%
%   Outputs:
%      qsiTot_flow = NCells x NFracs matrix of fractional volumetric 
%                    transport rate per meter width due to flow [m3/s/m]
%      ThetaCrit_i = NCells x NFracs matrix of critical shields stress in
%                    each cell, for each fraction (provided for use in
%                    calculation of bed slope effects on transverse
%                    transport rate).
%
%   References:
%      Meyer-Peter, E. & Müller, R., 1948. Formulas for bed-load transport.
%         In Proceedings of the 2nd Meeting of the International 
%         Association for Hydraulic Structures Research. pp. 39–64.
%      Wong, M. & Parker, G., 2006. Reanalysis and Correction of Bed-Load
%         Relation of Meyer-Peter and Müller Using Their Own Database.
%         Journal of Hydraulic Engineering, 132(11), pp.1159-1168.
%
%   See also: BEDLOAD, HIDINGFUNCTION, XCHANNELMODEL, WILCOCKCROWE.

ThetaCrit_i = Inputs.ST.ThetaCrit * HidExp;

Delta = (Inputs.Sed.Rho_S-Inputs.Hyd.Rho_W)/Inputs.Hyd.Rho_W;
            
qsiTot_flow = ...
    Cell.Fi * ...
    Inputs.ST.MPMcoef .* (ones(Cell.NCells,1)*Frac.Di_m) .* ...
    sqrt(Delta * Inputs.Hyd.g * (ones(Cell.NCells,1)*Frac.Di_m)) .* ...
    max(Cell.ShieldsStress_i - ThetaCrit_i,0).^Inputs.ST.MPMexponent;
end

