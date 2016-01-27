function [Delta_store] = StoreErosion(Inputs, Cell, Bank)
% Calculate rate of change in stored erosion volume over timestep [m3/m/s]
% Erosion stored until sufficient available to bank top cell to level of
% bank toe cell.
Delta_store = zeros(Cell.NCells,1);
for BankNo = 1:Bank.NBanks
    Top = Bank.Top(BankNo);
    Bottom = Bank.Bottom(BankNo);
    
    if (Cell.Z(Top) + ...
       (Cell.EroStore(Top) + Cell.Delta_bank(Top)) / Cell.Width(Top) / Inputs.Sed.Porosity) < ...
       Cell.Z(Bottom)
        % implement sudden erosion
        Delta_store(Top) = Cell.EroStore(Top) / Inputs.Time.dT;
    else
        % store erosion
        Delta_store(Top) = -Cell.Delta_bank(Top);
    end
    
end

end

