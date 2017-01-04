function [qsiTot_flow, ThetaCrit_i] = BedLoad(Inputs, Cell, Frac)
%BEDLOAD   Calculate fractional volumetric bedload transport due to flow
%BEDLOAD applies the user specified transport formula to calculate cell
%center transport rate. For simulations with multiple sediment fractions a
%fractional transport rate is calculated taking into account
%hiding/exposure using the user specified hiding function.
%
%   [qsiTot_flow, ThetaCrit_i] = BEDLOAD(Inputs, Cell, Frac)
%
%   References:
%      Meyer-Peter, E. & Müller, R., 1948. Formulas for bed-load transport.
%         In Proceedings of the 2nd Meeting of the International 
%         Association for Hydraulic Structures Research. pp. 39–64.
%      Wilcock, P.R. & Crowe, J.C., 2003. Surface-based transport model for
%         mixed-size sediment. Journal of Hydraulic Engineering, 129(2), 
%         pp.120–128.
%
%   See also: XCHANNELMODEL, HIDINGFUNCTION, MEYERPETERMULLER, WILCOCKCROWE

%% Apply hiding fuction to calculate fractional critical shear stress
HidExp = HidingFunction(Inputs.ST, Cell, Frac);

%% Apply sediment transport formula
switch Inputs.ST.Formula
    case 1
        % Meyer-Peter-Muller
        [qsiTot_flow, ThetaCrit_i] = MeyerPeterMuller(Inputs, Cell, Frac, HidExp);
    case 2
        % Wilcock-crowe formula (Wilcock and Crowe 2003)
        [qsiTot_flow, ThetaCrit_i] = WilcockCrowe(Inputs, Cell, Frac, HidExp);
end
Cell.qsiTot_flow(~Cell.Active) = 0;

end

