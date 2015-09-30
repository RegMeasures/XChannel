function [Q] = CalcQ(CellWidth,BedLevel,WetLastTimestep,Slope,Roughness,DryFlc,WL)
% Calculate cross-section flow for given geometry, WL, etc.

WetCells = (WL-BedLevel>=DryFlc)|(WetLastTimestep&(WL-BedLevel>=DryFlc/2));
H = WetCells.*(WL-BedLevel);
Q = H.^(5/3).*CellWidth*Slope^0.5/Roughness;
Q = sum(Q);
end

