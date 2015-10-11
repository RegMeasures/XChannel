function [qsiN_slope] = BedSlope(Inputs, Edge, Frac)
% Calc transverse fractional volumetric transport rate due to bed slope

% Note: at the moment, for the purposes of slope analysis we are
%       assuming unadjusted bedload transport vector is normal to
%       cross-section and that downstream slope << transverse slope
%       (i.e. for the bed slope calc secondary flow and downstream 
%       slope are neglected). 

switch Inputs.Opt.Slope.Formula
    case 0
        % No bed slope influence on transport
        qsiN_slope = zeros(Edge.NEdges, Frac.NFracs);
    case 1
        % General transverse bed slope formula in the form of Sekine & Parker 1992
        Beta = Inputs.Opt.Slope.BetaStar * (Edge.ThetaCrit_i./Edge.ShieldsStress_i).^Inputs.Opt.Slope.m;
        Beta(Edge.ShieldsStress_i == 0) = 0;
        qsiN_slope = - Edge.qsiTot_flow .* Beta .* (Edge.Slope * ones(1, Frac.NFracs));
    case 2
        % Talmon et al 1995
        
        fTheta = Inputs.Opt.Slope.A_sh * Edge.ShieldsStress_i.^Inputs.Opt.Slope.B_sh .* ...
                 (((ones(Edge.NEdges, 1) * Frac.Di_m)) ./ (Edge.H * ones(1, Frac.NFracs))).^Inputs.Opt.Slope.C_sh .* ...
                 ((ones(Edge.NEdges,1) * Frac.Di_m) ./ (Edge.Dg_m * ones(1,Frac.NFracs))).^Inputs.Opt.Slope.D_sh;
        qsiN_slope = - Edge.qsiTot_flow .* (1 ./ fTheta) .* (Edge.Slope * ones(1, Frac.NFracs));
        qsiN_slope(isnan(qsiN_slope)) = 0;
end

end

