function [qsiTot_flow,ThetaCrit_i] = MeyerPeterMuller(Inputs, Cell, Frac, HidExp)
% Calculate volumetric transport rate per unit width using Meyer-Peter-Muller formula

ThetaCrit_i = Inputs.ST.ThetaCrit * HidExp;

Delta = (Inputs.Sed.Rho_S-Inputs.Hyd.Rho_W)/Inputs.Hyd.Rho_W;
            
qsiTot_flow = Cell.Fi * ...
              Inputs.ST.MPMcoef .* (ones(Cell.NCells,1)*Frac.Di_m) .* ...
              sqrt(Delta * Inputs.Hyd.g * (ones(Cell.NCells,1)*Frac.Di_m)) .* ...
              max(Cell.ShieldsStress_i - ThetaCrit_i,0).^Inputs.ST.MPMexponent;
end

