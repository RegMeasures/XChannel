function [Q] = CalcQ(CellWidth,BedLevel,WetLastTimestep,HydInputs,WL)
% Calculate cross-section flow for given geometry, WL, etc.

WetCells = (WL-BedLevel>=HydInputs.DryFlc)|(WetLastTimestep&(WL-BedLevel>=HydInputs.DryFlc/2));
H = WetCells.*(WL-BedLevel);
switch HydInputs.Roughness
    case 1 % Manning roughness
        Chezy = (H.^(1/6)) / HydInputs.ManningN;
    case 2 % Colebrook-white
        Chezy = 18 * log10(12 * H / HydInputs.ks); % eqn 2.34b sedimentation engineering or 9.56 D3D
end
Chezy = max(Chezy, 0.001);
Q = H .* CellWidth .* Chezy .* sqrt(H * HydInputs.Slope);
Q = sum(Q);
end

