function [Delta_i_bank] = BankFlux(Options, Cell, Edge, Frac, dT, Bank)
% Calculate fractional volumetric sediment flux due to bank erosion

%qsiN_Bank = zeros(Edge.NEdges,1);
Delta_i_bank = zeros(Cell.NCells,Frac.NFracs);

for BankNo = 1:Bank.NBanks
    Top = Bank.Top(BankNo);
    Bottom = Bank.Bottom(BankNo);
    switch Options.Approach
        case 0
            % No bank erosion flux
            Flux = 0;
        case 1
            % Excess slope bank erosion flux
            CellSeperation = abs(Cell.N(Top) - Cell.N(Bottom));
            ExcessHeight = (Cell.Z(Top) - Cell.Z(Bottom)) - ...
                           Options.Repose * CellSeperation;
            Flux = Options.SlipRatio * ExcessHeight/2 * ...
                   sum(Cell.Width([Top,Bottom]))/2; % Total volumetric sed flux as a result of individual bank eroding [m3/m/dT]
    end
    Flux_i = Flux * Cell.Fi(Top,:) ./ dT;
    Delta_i_bank(Top,:) = Delta_i_bank(Top,:) - Flux_i;
    Delta_i_bank(Bottom,:) = Delta_i_bank(Bottom,:) + Flux_i;
end

end

