function [qsiTot_flow, ThetaCrit_i] = BedLoad(Inputs, Cell, Frac)
%BEDLOAD   Calculate fractional volumetric bedload transport due to flow
%BEDLOAD applies the user specified transport formula to calculate cell
%center transport rate. For simulations with multiple sediment fractions a
%fractional transport rate is calculated taking into account
%hiding/exposure using the user specified hiding function.
%
%   [qsiTot_flow, ThetaCrit_i] = BEDLOAD(Inputs, Cell, Frac)
%
%   Inputs:
%      Inputs = Struct of model inputs read in by ReadModelInputs.
%               Fields used by BEDLOAD are:
%         .ST.HidExp       = Choice of hiding/exposure formulation
%         .ST.Gamma        = Hiding function exponent
%         .ST.Formula      = Choice of bedload transport formula: 
%                            1-> Meyer-Peter-Muller (1948)
%                            2-> Wilcock and Crowe (2003)
%         .ST.ThetaCrit    = Critical shields stress (for MPM)
%         .ST.MPMcoef      = MPM coeffieint (for MPM)
%         .ST.MPMexponent  = MPM exponent (for MPM)
%         .Hyd.g           = Acceleration due to gravity [m/s2]
%         .Sed.Rho_S       = Sediment density [kg/m3]
%         .Hyd.Rho_W       = Water density [kg/m3]
%      Cell   = Struct of cell center properties including:
%         .NCells          = number of cells in cross-section
%         .Dg_m            = NCells x 1 matrix of geometric mean grain size
%                            of each cell [m]
%         .Fi              = NCells x NFrac matrix of proportion of active
%                            layer made up of each fraction, in each cell.
%         .ShieldsStress_i = NCells x NFracs matrix of shields stress
%                            affecting each fraction, in each cell.
%      Frac   = Struct of sediment fraction properties including:
%         .NFracs          = number of sediment fractions (integer)
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
%      Wilcock, P.R. & Crowe, J.C., 2003. Surface-based transport model for
%         mixed-size sediment. Journal of Hydraulic Engineering, 129(2), 
%         pp.120–128.
%
%   See also: XCHANNEL, HIDINGFUNCTION, MEYERPETERMULLER,
%   WILCOCKCROWE.

%% Apply hiding fuction to calculate fractional critical shear stress
HidExp = HidingFunction(Inputs.ST, Cell, Frac);

%% Apply sediment transport formula
switch Inputs.ST.Formula
    case 1
        % Meyer-Peter-Muller (1948)
        [qsiTot_flow, ThetaCrit_i] = MeyerPeterMuller(Inputs, Cell, Frac, HidExp);
    case 2
        % Wilcock-crowe formula (Wilcock & Crowe, 2003)
        [qsiTot_flow, ThetaCrit_i] = WilcockCrowe(Inputs, Cell, Frac, HidExp);
end
Cell.qsiTot_flow(~Cell.Active) = 0;

end

