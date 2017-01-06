function [Chezy, H, U, Tau_S, WetCells] = ...
    BasicHydraulics(BedLevel, WetLastTimestep, HydInputs, WL)
%BASICHYDRAULICS Calculate depth, Chezy roughness, velocity & shear stress
%
%   [Chezy, H, U, Tau_S, WetCells] = BASICHYDRAULICS(BedLevel, ...
%                                    WetLastTimestep, HydInputs, WL)
%
%   Inputs:
%      BedLevel  = NCells x 1 matrix of bed level of each cell [m]
%      WetLastTimestep = NCells x 1 boolean matrix of which cells were wet
%                        in the last timestep
%      HydInputs = Struct of user specified inputs relating to hydraulics
%                  as read in to Inputs.Hyd by ReadModelInputs. Fields of
%                  HydInputs used by BASICHYDRAULICS are:
%                    .Roughness = roughness formulation ID
%                    .ManningN or .ks = roughness coefficient [s/m^(1/3)]
%                                       or [m] respectively 
%                    .Slope     = streamwise bed slope [m/m]
%                    .DryFlc    = threshold depth for drying/flooding [m]
%      WL        = Water level [m]
%   Outputs:
%      Chezy     = NCells x 1 matrix of Chezy roughness coefficient
%      H         = NCells x 1 matrix of water depth [m]
%      U         = NCells x 1 matrix of depth averaged water velocity [m/s]
%      Tau_S     = NCells x 1 matrix of streamwise bed shear stress [N/m2]
%      WetCells  = NCells x 1 boolean matrix identifying wet cells
%
%   References:
%      Deltares, 2014. Delft3D-Flow User Manual v3.15.34158
%      García, M.H., 2008. Sediment Transport and Morphodynamics. In M. H.
%         García, ed. Sedimentation Engineering: Processes Measurements,
%         Modeling and Practice. American Society of Civil Engineers.
%   
%   See also: XCHANNEL, READMODELINPUTS, SETQ.

WetCells = (WL-BedLevel >= HydInputs.DryFlc) | ...
           (WetLastTimestep & (WL-BedLevel >= (HydInputs.DryFlc/2)));

H = WetCells .* (WL-BedLevel);

switch HydInputs.Roughness
    case 1 % Manning roughness
        Chezy = (H.^(1/6)) / HydInputs.ManningN;
    case 2 % Colebrook-white
        Chezy = 18 * log10(12 * H / HydInputs.ks);
        % eqn 2.34b sedimentation engineering or 9.56 D3D manual
end
Chezy = max(Chezy, 0.001);

U = Chezy .* sqrt(H * HydInputs.Slope);

Tau_S = HydInputs.Rho_W * HydInputs.g * H * HydInputs.Slope;

end

