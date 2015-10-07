function [Delta_i_bank] = BankFlux(Options, Cell, Edge, Frac, dT, Bank)
% Calculate fractional volumetric sediment flux rate due to bank erosion
% Delta_i_bank = NCells x NFracs matrix containing fractional flux rate
% into and out of each cell due to bank erosion.

%qsiN_Bank = zeros(Edge.NEdges,1);
Delta_i_bank = zeros(Cell.NCells,Frac.NFracs);

for BankNo = 1:Bank.NBanks
    Top = Bank.Top(BankNo);
    Bottom = Bank.Bottom(BankNo);
    switch Options.Approach
        case 0
            % No bank erosion flux
            FluxRate = 0;
        case 1
            % Excess slope bank erosion flux
            CellSeperation = abs(Cell.N(Top) - Cell.N(Bottom));
            ExcessHeight = (Cell.Z(Top) - Cell.Z(Bottom)) - ...
                           Options.Repose * CellSeperation;
            FluxRate = Options.SlipRatio * ExcessHeight/2 * ...
                       sum(Cell.Width([Top,Bottom]))/2 / dT; % Total volumetric sed flux as a result of individual bank eroding [m3/m/dT]
        case 2
            % Flux based on bank toe erosion rate (THETSD approach)
            ToeErosionRate = -sum(Cell.Delta_i_flow(Bottom,:) + Cell.Delta_i_slope(Bottom,:), 2);
            ToeErosionRate = max(ToeErosionRate,0);
            FluxRate = Options.ThetSD * ToeErosionRate;
        case 3
            % Flux based on bank toe transport rate
            ToeTransportRate = Cell.qsS_flow_kg(Bottom);
            FluxRate = Options.QsBeRatio * ToeTransportRate;
        case 4
            % Flux based on toe transport rate * slope
            Slope = (Cell.Z(Top) - Cell.Z(Bottom)) / ...
                    abs(Cell.N(Top) - Cell.N(Bottom));
            ToeTransportRate = Cell.qsS_flow_kg(Bottom);
            FluxRate = Options.BErodibility * ToeTransportRate * Slope;
    end
    FluxRate_i = FluxRate * Cell.Fi(Top,:);
    Delta_i_bank(Top,:) = Delta_i_bank(Top,:) - FluxRate_i; % MIGHT NEED TO CHANGE THIS BIT TO DEAL WITH SED MIXING WHEN TOP AND BOTTOM ARE NOT ADJACENT...?
    Delta_i_bank(Bottom,:) = Delta_i_bank(Bottom,:) + FluxRate_i;
end

end

