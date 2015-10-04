function [WL] = SetQ(CellWidth,BedLevel,WetLastTimestep,HydInputs)
% Find water level (normal depth) which generates the specified flow

% Use simple but robust bisection method
WL = bisection(@Qerr,max(BedLevel),min(BedLevel),HydInputs.QTol);
    function y = Qerr(x)
        y=CalcQ(CellWidth,BedLevel,WetLastTimestep,HydInputs.Slope,...
            HydInputs.Roughness,HydInputs.DryFlc,x)-HydInputs.Flow;
    end

% Use inbuilt matlab optimisation (requires optimization toolbox)
% WL = fsolve(@Qerr,1);
%     function y = Qerr(x)
%         y=Flow-CalcQ(CellWidth,BedLevel,Slope,Roughness,DRYFLC,x);
%     end

end

