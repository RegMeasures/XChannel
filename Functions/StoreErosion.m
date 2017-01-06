function [Delta_store] = StoreErosion(Inputs, Cell, Bank, dT)
%STOREEROSION   Rate of change in stored erosion volume over timestep
%STOREEROSION implements the Nicholas (2013) approach of maintaining bank
%height by storing bank erosion until sufficient to lower bank top cell to
%height of bank toe cell.
%
%   Inputs:
%      Cell       = Struct of cell center properties initialised by
%                   InitialiseVariables and set in earlier steps of 
%                   XChannel
%      Bank       = Struct of bank properties created by BankStencil
%      dT         = model timestep [s]
%   
%   Outputs:
%      Delta_store = NCells x 1 matrix of rate of change of 
%                    'stored bank erosion' volume in each cell [m3/m/s].
%                    Delta_store is positive when erosion is being 'stored'
%                    and negative when it is being 'released'.
%
%   References:
%      Nicholas, A.P., 2013. Modelling the continuum of river channel
%         patterns. Earth Surface Processes and Landforms.
%   
%   See also: XCHANNEL, IDENTIFYBANKS, BANKSTENCIL, TRIGGERBANKS, BANKFLUX,
%   INITIALISEVARIABLES, READMODELINPUTS

Delta_store = zeros(Cell.NCells,1);

% Loop through banks
for BankNo = 1:Bank.NBanks 
    Top = Bank.Top(BankNo);
    Bottom = Bank.Bottom(BankNo);
    
    % If stored bank erosion volume combined with erosion this timestep is 
    % sufficient to reduce bank height to toe elevation
    % note: EroStore and Delta_bank(Top) variables are both negative
    if (Cell.Z(Top) + ...
        (Cell.EroStore(Top) + Cell.Delta_bank(Top)*dT) / ...
        (Cell.Width(Top) * (1-Inputs.Sed.Porosity))) < ...
       Cell.Z(Bottom)
        % then implement sudden erosion
        Delta_store(Top) = Cell.EroStore(Top) / Inputs.Time.dT;
    else
        % otherwise store erosion
        Delta_store(Top) = -Cell.Delta_bank(Top);
    end
    
end

end

