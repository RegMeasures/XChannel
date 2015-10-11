function [qsiTot_flow, ThetaCrit_i] = BedLoad(Inputs, Cell, Frac)
% Calculate fractional volumetric bedload transport due to flow (in cell centres)

%% Apply hiding and exposure fuction to calculate fractional critical shear stress
HidExp = HidingFunction(Inputs.Opt.ST, Cell, Frac, Cell.Dg_m);

%% Apply sediment transport formula
switch Inputs.Opt.ST.Formula
    case 1
        % Meyer-Peter-Muller
        [qsiTot_flow, ThetaCrit_i] = MeyerPeterMuller(Inputs, Cell, Frac, HidExp);
    case 2
        % Wilcock-crowe formula
        [qsiTot_flow, ThetaCrit_i] = WilcockCrowe(Inputs, Cell, Frac, HidExp);
end
Cell.qsiTot_flow(~Cell.Active) = 0;

end

