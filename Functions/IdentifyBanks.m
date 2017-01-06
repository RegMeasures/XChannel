function [IsBank] = IdentifyBanks(Options, Cell, Edge)
%IDENTIFYBANKS   ID banks using selected bank identification routine
%
%   [IsBank] = IDENTIFYBANKS(Options, Cell, Edge)
%
%   Inputs:
%      Options = Bank identification calculation options as read in to 
%                Inputs.Bank.Stencil struct by ReadModelInputs. Fields are:
%         .Approach = Choice of bank identification approach:
%                     Approach = 0 -> everywhere is a bank
%                     Approach = 1 -> wet/dry bank ID
%                     Approach = 2 -> transport/no-transport bank ID
%                     Approach = 3 -> bank height bank ID
%                     Approach = 4 -> bank slope bank ID
%         .BHeight  = Height threshold [m] (req'd for Approach = 3)
%         .BSlope   = Slope threshold [m/m] (req'd for Approach = 4)
%      Cell    = Struct of cell center properties initialised by
%                InitialiseVariables and set in earlier steps of 
%                XChannel
%      Edge    = Struct of cell edge properties initialised by
%                InitialiseVariables and set in earlier steps of 
%                XChannel
%
%   Outputs:
%      IsBank = NEdges x 1 matrix of bank locations: Left banks = +1, Right
%               banks = -1, Non-bank edges = 0.
%
%   See also: XCHANNEL, BANKSTENCIL, TRIGGERBANKS, BANKFLUX,
%   INITIALISEVARIABLES, READMODELINPUTS.

switch Options.Approach
    case 0 
        % Everywhere is a bank (no bank ID routine)
        % note: This is not as trivial as it seems because we still
        % identify bank direction.
        BankHeight = Cell.Z(2:end) - Cell.Z(1:end-1);
        IsBank = (BankHeight <= 0) - (BankHeight > 0);
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
        IsBank = (BankHeight <= -Options.BHeight) - ...
                 (BankHeight >= Options.BHeight);
    case 4
        % Bank slope bank identification approach
        IsBank = (Edge.Slope(2:end-1) <= -Options.BSlope) - ...
                 (Edge.Slope(2:end-1) >= Options.BSlope);
end

% Pad with zeros to indicate no bank erosion at cross-section ends
IsBank = [0; IsBank; 0];

end

