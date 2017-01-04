function [Tau_N] = SecondaryFlow(HydInputs, Cell)
%SECONDARYFLOW   Transverse shear stress due to secondary (spiral) flow
%Calculate transverse shear stess due to secondary flow (also known as
%spiral motion) induced by channel curvature. SECONDARYFLOW assumes
%fully developed secondary flow induced by constant radius streamline
%curvature.
%
%   [Tau_N] = SecondaryFlow(HydInputs, Cell)
%
%   Inputs:
%
%   Outputs:
%
%   References:
%      Deltares, 2014. Delft3D-Flow User Manual Version 3.15.34158.
%      Kalkwijk, J.P.T. & Booij, R., 1986. Adaptation of secondary flow in 
%         nearly-horizontal flow. Journal of Hydraulic Research, 24(1), 
%         pp.19–37.
%
%   See also: XCHANNELMODEL
%% Initialise variables
AlphaSpiral = NaN(Cell.NCells,1);
SpiralIntensity = zeros(Cell.NCells,1);
Tau_N = zeros(Cell.NCells,1);

%% Calculate secondary flow according to Kalkwijk & Booij 1986
AlphaSpiral(Cell.Wet) = min(sqrt(HydInputs.g) ./ ...
                            (HydInputs.Kappa * Cell.Chezy(Cell.Wet)), 0.5); 
% Equation 9.156 in Delft3D-FLOW User Manual or Eq8 Kalkwijk & Booij 1986

SpiralIntensity(Cell.Wet) = Cell.H(Cell.Wet) .* Cell.U(Cell.Wet) / ...
                            HydInputs.Radius; 
% Eq 9.160 in Delft3D-FLOW User Manual (With I=I_be (intensity due to bend)
% as I_ce (coriolis) is negligable)

Tau_N(Cell.Wet) = -2 * HydInputs.ESpiral * ...
                  HydInputs.Rho_W * ...
                  AlphaSpiral(Cell.Wet).^2 .* ...
                  (1 - AlphaSpiral(Cell.Wet) / 2) .* ...
                  Cell.U(Cell.Wet) .* ...
                  SpiralIntensity(Cell.Wet); 
% Eq9.145 in Delft3D-FLOW User Manual or Eq26 Kalkwijk & Booij 1986

end

