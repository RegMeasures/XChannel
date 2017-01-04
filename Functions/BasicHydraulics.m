function [Chezy, H, U, Tau_S, WetCells] = BasicHydraulics(BedLevel, WetLastTimestep, HydInputs, WL)
%BASICHYDRAULICS Calculate depth and Chezy roughness
%
%   [Chezy, H, U, Tau_S, WetCells] = BASICHYDRAULICS(BedLevel, ...
%                                    WetLastTimestep, HydInputs, WL)

WetCells = (WL-BedLevel >= HydInputs.DryFlc) | ...
           (WetLastTimestep & (WL-BedLevel >= (HydInputs.DryFlc/2)));

H = WetCells .* (WL-BedLevel);

switch HydInputs.Roughness
    case 1 % Manning roughness
        Chezy = (H.^(1/6)) / HydInputs.ManningN;
    case 2 % Colebrook-white
        Chezy = 18 * log10(12 * H / HydInputs.ks); % eqn 2.34b sedimentation engineering or 9.56 D3D
end
Chezy = max(Chezy, 0.001);

U = Chezy .* sqrt(H * HydInputs.Slope);

Tau_S = HydInputs.Rho_W * HydInputs.g * H * HydInputs.Slope;

end

