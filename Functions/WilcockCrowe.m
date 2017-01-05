function [qsiTot_flow,ThetaCrit_i] = WilcockCrowe(Inputs, Cell, Frac, ...
                                                  HidExp)
%WILCOCKCROWE   Bedload transport rate using Wilcock and Crowe (2003)
%Fractional volumetric transport rate per unit width calculated using the
%Wilcock and Crowe (2003) formulation. Note that this function uses
%externally calculated hiding/exposure corrections but for consistency
%these should be calculated using the Wilcock & Crowe approach (i.e. 
%HidExp = 2 in HidingFunction).
%
%   [qsiTot_flow, ThetaCrit_i] = WILCOCKCROWE(Inputs, Cell, Frac, HidExp)
%
%   Inputs:
%      Inputs = Struct of model inputs read in by ReadModelInputs.
%               Fields used by MEYERPETERMULLER are:
%         .Hyd.g     = Acceleration due to gravity [m/s2]
%         .Sed.Rho_S = Sediment density [kg/m3]
%         .Hyd.Rho_W = Water density [kg/m3]
%      Cell   = Struct of cell center properties including:
%         .NCells    = number of cells in cross-section
%         .Fi        = NCells x NFrac matrix of proportion of active
%                      layer made up of each fraction, in each cell.
%         .Tau_Tot   = NCells x 1 matrix of total shear stress
%                      affecting each each cell.
%      Frac   = Struct of sediment fraction properties including:
%         .Di_m      = 1 x NFracs matrix of fraction sediment size[m]
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
%      Wilcock, P.R. & Crowe, J.C., 2003. Surface-based transport model for
%         mixed-size sediment. Journal of Hydraulic Engineering, 129(2), 
%         pp.120–128.
%
%   See also: BEDLOAD, HIDINGFUNCTION, XCHANNEL, MEYERPETERMULLER.

Fs = sum(Cell.Fi(:,Frac.SandFrac),2);
Delta = (Inputs.Sed.Rho_S-Inputs.Hyd.Rho_W)/Inputs.Hyd.Rho_W;


TauStar_rm = 0.021 + 0.015 * exp(-20 * Fs);
Tau_rm = TauStar_rm * (Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * Inputs.Hyd.g .* Cell.Dg_m;

DiDg = (ones(Cell.NCells,1) * Frac.Di_m) ./ ...
       (Cell.Dg_m * ones(1,Frac.NFracs));
% Note: from HidingFunction with (HidExp = 2):
%    b = 0.67 ./ (1 + exp(1.5 - DiDg));
%    HidExp = DiDg.^b;
Tau_ri = HidExp .* (Tau_rm*ones(1,Frac.NFracs));
Phi = (Cell.Tau_Tot * ones(1,Frac.NFracs)) ./ Tau_ri;

WiStar = zeros(Cell.NCells, Frac.NFracs);
WiStar(Phi<1.35) = 0.002 * Phi(Phi<1.35).^7.5;
WiStar(Phi>=1.35) = 14 * (1 - 0.894 ./ Phi(Phi>=1.35).^0.5).^4.5;

qsiTot_flow = WiStar .* Cell.Fi .* ...
              (Cell.UStar*ones(1,Frac.NFracs)).^3 ./ ...
              (Delta * Inputs.Hyd.g);

Cell.ShieldsStress_i = (Cell.Tau_Tot * ones(1,Frac.NFracs)) ./ ...
                       ((Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * ...
                        Inputs.Hyd.g * ...
                        (ones(Cell.NCells,1) * Frac.Di_m));

ThetaCrit_i = Tau_ri ./ ((Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * ...
                         Inputs.Hyd.g * ...
                         (ones(Cell.NCells,1) * Frac.Di_m)); 
% Note: reference shear stress used instead of critical not sure if this is
% significant but it is only used for lateral slope effect calculations?
end

