function [qsiTot_flow] = WilcockCrowe(Rho_S, Rho_W, G, NCells, NFracs, D50i_m, SandFrac, Fi, Dg_m, Tau_Tot, UStar)
% Calculate volumetric transport rate per unit width using Wilcock-Crowe formula

Fs = sum(Fi(:,SandFrac),2);
Delta = (Rho_S-Rho_W)/Rho_W;

b = 0.67 ./ (1 + exp(1.5 - (1 ./ Dg_m) * D50i_m));
TauStar_rm = 0.021 + 0.015 * exp(-20 * Fs);
Tau_rm = TauStar_rm * (Rho_S - Rho_W) * G .* Dg_m;
Tau_ri = ((1 ./ Dg_m) * D50i_m).^b .* (Tau_rm*ones(1,NFracs));
Phi = (Tau_Tot*ones(1,NFracs)) ./ Tau_ri;

WiStar = zeros(NCells, NFracs);
WiStar(Phi<1.35) = 0.002 * Phi(Phi<1.35).^7.5;
WiStar(Phi>=1.35) = 14 * (1 - 0.894 ./ Phi(Phi>=1.35).^0.5).^4.5;

qsiTot_flow = WiStar .* Fi .* (UStar*ones(1,NFracs)).^3 ./ (Delta * G);

end

