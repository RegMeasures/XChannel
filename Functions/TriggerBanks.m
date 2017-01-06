function [Active] = TriggerBanks(Options, Cell, Bank)
%TRIGGERBANKS   Identify actively eroding banks
%Identify which of the banks are actively eroding using a variety of
%optional criteria.
%
%   [Active] = TRIGGERBANKS(Options, Cell, Bank)
%
%   Inputs:
%      Options = Bank trigger options as read in to Inputs.Bank.Trigger
%                struct by ReadModelInputs. Fields are:
%         .BTrigger = Type of bank trigger selected:
%                     BTrigger = 0 -> every bank is active
%                     BTrigger = 1 -> bank active if > threshold height
%                     BTrigger = 2 -> bank active if toe degrading
%                     BTrigger = 3 -> bank active if > threshold slope
%         .BTHeight = Height threshold [m] (for BTrigger = 1)
%         .BTSlope  = Slope threshold [m/m] (for BTrigger = 3)
%      Cell    = Struct of cell center properties initialised by
%                InitialiseVariables and set in earlier steps of 
%                XChannel
%      Bank    = Struct of bank properties initialised by
%                InitialiseVariables and set in earlier steps of 
%                XChannel
%   Outputs:
%      Active = NBanks x 1 boolean matrix indicating which banks are active
%
%   See also: XCHANNEL, IDENTIFYBANKS, BANKSTENCIL, BANKFLUX,
%   INITIALISEVARIABLES, READMODELINPUTS.

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
        Active = (Cell.Delta_flow(Bank.Bottom) + ...
                  Cell.Delta_slope(Bank.Bottom)) < 0;
    case 3
        % Banks active if bank slope exceeds threshold
        BankSlope = (Cell.Z(Bank.Top) - Cell.Z(Bank.Bottom)) / ...
                    abs(Cell.N(Bank.Top) - Cell.N(Bank.Bottom));
        Active = BankSlope > Options.BTSlope;
end
end