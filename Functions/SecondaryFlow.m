function [Tau_N] = SecondaryFlow(HydInputs, Cell)
%SECONDARYFLOW   Transverse shear stress due to secondary (spiral) flow
%Calculate transverse (normal) shear stess due to secondary flow (also 
%known as spiral motion) induced by channel curvature. SECONDARYFLOW 
%assumes fully developed secondary flow induced by constant radius 
%streamline curvature.
%
%   [Tau_N] = SecondaryFlow(HydInputs, Cell)
%
%   Inputs:
%      HydInputs = Struct of user specified inputs relating to hydraulics 
%                  as read in to Inputs.Hyd by ReadModelInputs. Fields of
%                  HydInputs used by SECONDARYFLOW are:
%                     .g       = acceleration due to gravity [m/s2]
%                     .Kappa   = Von Karman's constant
%                     .Chezy   = Chézy friction coefficient [m^0.5/s]
%                     .Radius  = specified bend radius [m]
%                     .ESpiral = user specified coefficient for spiral
%                                motion
%                     .Rho_W   = water density [kg/m3]
%                     .ConstBend = Option to use constant bend radius [1/0]
%      Cell      = Struct of cell center properties initialised by
%                  InitialiseVariables and set in earlier steps of 
%                  XChannel. Fields of Cell used by SECONDARYFLOW are:
%                     .NCells = number of cells across cross-section
%                     .Wet    = mask indicating wet cells
%                     .U      = depth averaged velocity in each cell [m/s]
%                     .H      = water depth in each cell [m]
%                     .N      = cell center location, cross-channel dir [m]
%   
%   Outputs:
%      Tau_N     = Transverse shear stress due to secondary currents [N/m2]
%
%   References:
%      Deltares, 2014. Delft3D-Flow User Manual Version 3.15.34158.
%      Kalkwijk, J.P.T. & Booij, R., 1986. Adaptation of secondary flow in 
%         nearly-horizontal flow. Journal of Hydraulic Research, 24(1), 
%         pp.19–37.
%
%   See also: XCHANNEL, READMODELINPUTS.

%% Initialise variables
AlphaSpiral = NaN(Cell.NCells,1);
SpiralIntensity = zeros(Cell.NCells,1);
Tau_N = zeros(Cell.NCells,1);

%% Calculate secondary flow according to Kalkwijk & Booij 1986
AlphaSpiral(Cell.Wet) = min(sqrt(HydInputs.g) ./ ...
                            (HydInputs.Kappa * Cell.Chezy(Cell.Wet)), 0.5); 
% Equation 9.156 in Delft3D-FLOW User Manual or Eq8 Kalkwijk & Booij 1986

if HydInputs.ConstBend == 1
    % Costant bend radius
    SpiralIntensity(Cell.Wet) = Cell.H(Cell.Wet) .* Cell.U(Cell.Wet) / ...
                                HydInputs.Radius; 
else
    % Bend radius varies across section
    SpiralIntensity(Cell.Wet) = Cell.H(Cell.Wet) .* Cell.U(Cell.Wet) ./ ...
                                (Cell.N(Cell.Wet) + HydInputs.Radius); 
end
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

