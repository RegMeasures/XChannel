function [IsBank] = IdentifyBanks(Options, Cell, Edge)
% Identify bank using selected bank identification routine

switch Options.Approach
    case 0 
        % Everywhere is a bank (no bank ID routine)
        IsBank = ones(Cell.NCells-1, 1);
    case 1 
        % Wet/dry bank identification
        IsBank = Cell.Wet(2:end) - Cell.Wet(1:end-1);
    case 2
        % Transporting/non-transporting bank identification
        Transporting = Cell.qsS_flow_kg > 0;
        IsBank = Transporting(2:end) - Transporting(1:end-1);
    case 3
        % Bank height bank identification approach
        BankHeight = Cell.Z(2:end) - Cell.Z(1:end-1);
        IsBank = (BankHeight >= Options.BHeight) - (BankHeight <= -Options.BHeight);
    case 4
        % Bank slope bank identification approach
        IsBank = (Edge.Slope(2:end-1) <= -Options.BSlope) - (Edge.Slope(2:end-1) >= Options.BSlope);
end

% Pad with zeros to indicate no bank erosion at section ends
        IsBank = [0; IsBank; 0];

end

