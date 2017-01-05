function [HidExp] = HidingFunction(Options, Cell, Frac)
%HIDINGFUNCTION   Hiding/exposure function for grain size mixtures
%HIDINGFUNCTION returns a critical shear stress multiplier for each grain
%size fraction in a mixture, for each computational cell.
%
%   [HidExp] = HidingFunction(Options, Cell, Frac)
%   
%   Inputs:
%      Options = Struct of sediment transport options as read in to
%                Inputs.ST by ReadModelInputs. Includes:
%         .HidExp = Choice of hiding/exposure formulation
%         .Gamma  = Hiding function exponent for Parker Klingeman & McLean
%                   (1982). Not required for other formulations.
%      Cell    = Struct of cell center properties including:
%         .NCells = number of cells in cross-section
%         .Dg_m   = NCells x 1 matrix of geometric mean grain size of each 
%                   cell [m]
%      Frac    = Struct of sediment fraction properties including:
%         .NFracs = number of sediment fractions (integer)
%         .Di_m   = 1 x NFracs matrix of fraction sediment size [m]
%   
%   Outputs:
%      HidExp  = NCells x NFracs matrix of critical shear stress 
%                multipliers for each grain size fraction in a mixture, for
%                each computational cell. Theta_crit_i/Theta_crit = HidExp.
%   
%   References:
%      Ashida, K. & Michiue, M., 1972. Study on hydraulic resistance and
%         bedload transport rate in alluvial streams. Transactions, Japan
%         Society of Civil Engineering (In Japanese), 206 pp.59-69.
%      Parker, G., Klingeman, P. & McLean, D., 1982. Bedload and Size
%         Distribution in Paved Gravel-Bed Streams. Journal of the 
%         Hydraulics Division ASCE, 108(4) pp.544-571.
%      Wilcock, P.R. & Crowe, J.C., 2003. Surface-based transport model for
%         mixed-size sediment. Journal of Hydraulic Engineering, 129(2), 
%         pp.120–128.
%   
%   See also: XCHANNELMODEL, BEDLOAD, READMODELINPUTS.

switch Options.HidExp
    case 0
        % No hiding or exposure
        HidExp = ones(Cell.NCells, Frac.NFracs);
    case 1
        % Ashida and Michiue (1972)
        DiDg = (ones(Cell.NCells,1) * Frac.Di_m) ./ ...
               (Cell.Dg_m * ones(1,Frac.NFracs));
        HidExp = 0.8429 ./ DiDg;
        HidExp(DiDg >= 0.38889) = ...
            (log10(19) ./ (log10(19) + log10(DiDg(DiDg >= 0.38889)))).^2;
    case 2
        % Wilcock and Crowe (2003)
        DiDg = (ones(Cell.NCells,1) * Frac.Di_m) ./ ...
               (Cell.Dg_m * ones(1,Frac.NFracs));
        b = 0.67 ./ (1 + exp(1.5 - DiDg));
        HidExp = DiDg.^b;
    case 3
        % Parker Klingeman & McLean (1982)
        HidExp =  ((Cell.Dg_m * ones(1,Frac.NFracs)) ./ ...
                   (ones(Cell.NCells,1) * Frac.Di_m) ).^Options.Gamma;
end

end

