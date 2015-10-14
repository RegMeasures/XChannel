function [HidExp] = HidingFunction(Options, Cell, Frac)
% Hiding function for calculation of critical shields stress in grain size mixtures

switch Options.HidExp
    case 0
        % No hiding or exposure
        HidExp = ones(Cell.NCells,Frac.NFracs);
    case 1
        % Ashida and Michiue
        DiDg = (ones(Cell.NCells,1)*Frac.Di_m) ./ (Cell.Dg_m*ones(1,Frac.NFracs));
        HidExp = 0.8429 ./ DiDg;
        HidExp(DiDg >= 0.38889) = (log10(19) ./ (log10(19) + log10(DiDg(DiDg >= 0.38889)))).^2;
    case 2
        % Wilcock and Crowe
        DiDg = (ones(Cell.NCells,1)*Frac.Di_m) ./ (Cell.Dg_m*ones(1,Frac.NFracs));
        b = 0.67 ./ (1 + exp(1.5 - DiDg));
        HidExp = DiDg.^(b-1);
    case 3
        % Parker Klingeman & McLean
        HidExp =  ((Cell.Dg_m*ones(1,Frac.NFracs)) ./ (ones(Cell.NCells,1)*Frac.Di_m)).^Options.Gamma;
end

end

