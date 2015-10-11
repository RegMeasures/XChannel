function [HidExp] = HidingFunction(Options, Cell, Frac, Dg_m)
% Hiding function for calculation of critical shields stress in grain size mixtures

switch Options.HidExp
    case 0
        % No hiding or exposure
        HidExp = ones(Cell.NCells,Frac.NFracs);
    case 1
        % Ashida and Michiue
        DiDm = (ones(Cell.NCells,1)*Frac.Di_m) ./ (Dg_m*ones(1,Frac.NFracs));
        HidExp = 0.8429 ./ DiDm;
        HidExp(DiDm >= 0.38889) = (log10(19) ./ (log10(19) + log10(DiDm(DiDm >= 0.38889)))).^2;
    case 2
        % Wilcock and Crowe
        DiDm = (ones(Cell.NCells,1)*Frac.Di_m) ./ (Dg_m*ones(1,Frac.NFracs));
        b = 0.67 ./ (1 + exp(1.5 - DiDm));
        HidExp = DiDm.^b;
end

end

