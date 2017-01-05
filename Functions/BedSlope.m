function [qsiN_slope] = BedSlope(Slope, Edge, Frac)
%BEDSLOPE   Transverse transport rate due to bed slope effect
% BEDSLOPE calculates the fractional volumetric transverse transport rate
% due to transverse bed slope effects on fluvial transport.
%
%   [qsiN_slope] = BEDSLOPE(Slope, Edge, Frac)
%
%   Notes: 
%      For the purposes of slope analysis this function assumes the 
%      unadjusted bedload transport vector is normal to the cross-section 
%      and that downstream slope << transverse slope
%      (i.e. secondary flow and downstream slope are neglected). 
%
%   References:
%      Sekine, M. & Parker, G., 1992. Bed-load transport on transverse 
%         slope I. Journal of Hydraulic Engineering, 118(4), pp.513-535.
%      Talmon, A.M., Struiksma, N. & Van Mierlo, M.C.L.., 1995. Laboratory 
%         measurements of the direction of sediment transport on transverse
%         alluvial-bed slopes. Journal of Hydraulic Research, 33(4), 
%         pp.495-517.
%
%   See also: XCHANNEL, BEDLOAD

switch Slope.Formula
    case 0
        % No bed slope influence on transport
        qsiN_slope = zeros(Edge.NEdges, Frac.NFracs);
    case 1
        % General bed slope formula in the form of Sekine & Parker (1992)
        Beta = Slope.BetaStar * ...
               (Edge.ThetaCrit_i./Edge.ShieldsStress_i).^Slope.m;
        Beta(Edge.ShieldsStress_i == 0) = 0;
        qsiN_slope = - Edge.qsiTot_flow .* Beta .* ...
                     (Edge.Slope * ones(1, Frac.NFracs));
    case 2
        % Talmon et al (1995) transverse bed slope formula
        fTheta = Slope.A_sh * ...
                 Edge.ShieldsStress_i.^Slope.B_sh .* ...
                 (((ones(Edge.NEdges, 1) * Frac.Di_m)) ./ ...
                  (Edge.H * ones(1, Frac.NFracs))).^Slope.C_sh .* ...
                 ((ones(Edge.NEdges,1) * Frac.Di_m) ./ ...
                  (Edge.Dg_m * ones(1,Frac.NFracs))).^Slope.D_sh;
        qsiN_slope = - Edge.qsiTot_flow .* (1 ./ fTheta) .* ...
                     (Edge.Slope * ones(1, Frac.NFracs));
        qsiN_slope(isnan(qsiN_slope)) = 0;
end

end

