function [Active] = TriggerBanks(Options, Cell, Bank)
% Identify active banks using appropriate bank trigger
%

% Apply active bank trigger
switch Options.BTrigger
    case 0
        % Every bank is active
        Active = true(Bank.NBanks,1);
    case 1
        % Banks active if above threshold height
        BankHeight = Cell.Z(Bank.Top) - Cell.Z(Bank.Bottom);
        Active = BankHeight > Options.BTHeight;
    case 2
        % Banks active if bank toe is degrading due to slope and flow
        % (without effect of bank erosion)
        Active = (Cell.Delta_flow(Bank.Bottom) + Cell.Delta_slope(Bank.Bottom)) < 0;
    case 3
        % Banks active if bank slope exceeds threshold
        BankSlope = (Cell.Z(Bank.Top) - Cell.Z(Bank.Bottom)) / ...
                    abs(Cell.N(Bank.Top) - Cell.N(Bank.Bottom));
        Active = BankSlope > Options.BTSlope;
end
end