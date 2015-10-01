function [EdgeVals] = Centre2Edge(CentreVals, C2E_Weights)
% Convert cell centre parameter to sell edge parameter using C2E_Weights to
% apply either upwinding or central approaches

NCols = size(CentreVals,2);
C2E_Weights = C2E_Weights*ones(1, NCols);
EdgeVals = [zeros(1,NCols);
            CentreVals(1:end-1,:) .* C2E_Weights + CentreVals(2:end,:) .* (1-C2E_Weights); 
            zeros(1,NCols)];

end

