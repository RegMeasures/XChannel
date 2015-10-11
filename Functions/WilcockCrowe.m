function [qsiTot_flow,ThetaCrit_i] = WilcockCrowe(Inputs, Cell, Frac, HidExp)
% Calculate volumetric transport rate per unit width using Wilcock-Crowe formula

Fs = sum(Cell.Fi(:,Frac.SandFrac),2);
Delta = (Inputs.Sed.Rho_S-Inputs.Hyd.Rho_W)/Inputs.Hyd.Rho_W;

%b = 0.67 ./ (1 + exp(1.5 - (1 ./ Cell.Dg_m) * Frac.Di_m));
TauStar_rm = 0.021 + 0.015 * exp(-20 * Fs);
Tau_rm = TauStar_rm * (Inputs.Sed.Rho_S - Inputs.Hyd.Rho_W) * Inputs.Hyd.g .* Cell.Dg_m;
%Tau_ri = ((1 ./ Cell.Dg_m) * Frac.Di_m).^b .* (Tau_rm*ones(1,Frac.NFracs));
Tau_ri = HidExp .* (Tau_rm*ones(1,Frac.NFracs));
Phi = (Cell.Tau_Tot*ones(1,Frac.NFracs)) ./ Tau_ri;

WiStar = zeros(Cell.NCells, Frac.NFracs);
WiStar(Phi<1.35) = 0.002 * Phi(Phi<1.35).^7.5;
WiStar(Phi>=1.35) = 14 * (1 - 0.894 ./ Phi(Phi>=1.35).^0.5).^4.5;

qsiTot_flow = WiStar .* Cell.Fi .* (Cell.UStar*ones(1,Frac.NFracs)).^3 ./ (Delta * Inputs.Hyd.g);

ThetaCrit_i = Tau_ri; % for use in other places... reference shear stress used instead of critical.
end

