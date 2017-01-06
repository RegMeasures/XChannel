function [WL] = SetQ(CellWidth,BedLevel,WetLastTimestep,HydInputs)
%SETQ   Find water level (normal depth) which generates the specified flow
%Optimises WL in a cross-section so that total flow through the
%cross-section equals desired flow through the cross-section (for the
%specified geometry, roughness, slope and other options). A simple
%bisection algorithm is used to identify the WL which returns the desired
%flow.
%
%   [WL] = SETQ(CellWidth,BedLevel,WetLastTimestep,HydInputs)
%
%   Inputs:
%      CellWidth       = Width of each cell in cross-section [m]
%      BedLevel        = Bed elevation of each cell in cross-section [m]
%      WetLastTimestep = Mask of cells which were wet at the last timestep
%      HydInputs       = Struct of user specified inputs relating to 
%                        hydraulics as read in to Inputs.Hyd by 
%                        ReadModelInputs. Fields of HydInputs used by 
%                        SETQ are:
%                           .Flow      = desired flow [m3/s]
%                           .Roughness = roughness formulation ID
%                           .ManningN or .ks = roughness coefficient 
%                                        [s/m^(1/3)] or [m] respectively
%                           .Slope     = streamwise bed slope [m/m]
%                           .DryFlc    = threshold depth for drying/
%                                        flooding [m]
%                           .QTol      = acceptable error for flow [m3/s]
%                           .ItMax     = maximum number of iterations
%                                        allowed to resolve flow
%
%   Outputs:
%      WL              = Water level which provides the desired flow (m)
%
%   See also: XCHANNEL, BISECTION.

%% Use bisection method (simple & robust) to find WL which gives desired Q
WL = bisection(@Qerr, min(BedLevel), max(BedLevel),...
               HydInputs.QTol, HydInputs.ItMax);
    function y = Qerr(x)
        y = CalcQ(CellWidth,BedLevel,WetLastTimestep,HydInputs,x)...
            - HydInputs.Flow;
    end

%% An alternative to bisection is to use inbuilt matlab optimisation 
% (discarded as this requires optimization toolbox)

% WL = fsolve(@Qerr,1);
%     function y = Qerr(x)
%         y=Flow-CalcQ(CellWidth,BedLevel,Slope,Roughness,DRYFLC,x);
%     end

end

function [Q] = CalcQ(CellWidth, BedLevel, WetLastTimestep, HydInputs, WL)
%CALCQ   Calculate cross-section flow for given water level

[Chezy, H, ~, ~, ~] = BasicHydraulics(BedLevel, WetLastTimestep, ...
                                      HydInputs, WL);

Q = H .* CellWidth .* Chezy .* sqrt(H * HydInputs.Slope);
Q = sum(Q);
end