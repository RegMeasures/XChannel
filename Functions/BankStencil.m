function [Bank] = BankStencil(Options, Cell, Edge)
%BANKSTENCIL   Apply bank stencil to identify bank top and bottom (toe)
%For the purposes of bank erosion implementation in 2D models it is not 
%always appropriate to have bank top and bottom (toe) cells located
%adjacent to each other. BANKSTENCIL identifies the top and bottom cells
%associated with each identified bank. The identified cells are then used
%for flux calculation and bed level updating.
%
%   [Bank] = BANKSTENCIL(Options, Cell, Edge)
%
%   Inputs:
%      Options = Bank flux calculation options as read in to 
%                Inputs.Bank.Stencil struct by ReadModelInputs by
%      Cell    = Struct of cell center properties initialised by
%                InitialiseVariables and set in earlier steps of 
%                XChannelModel
%      Edge    = Struct of cell edge properties initialised by
%                InitialiseVariables and set in earlier steps of 
%                XChannelModel
%   Outputs:
%      Bank = Struct containing details of bank locations and top/bottom:
%         Bank.NBanks = Number of active banks currently in model
%         Bank.Top    = vector (size [NBanks,1]) listing cellNo for the top
%                       of each bank
%         Bank.Bottom = vector (size [NBanks,1]) listing cellNo for the 
%                       bottom of each bank
%
%   See also: XCHANNELMODEL, IDENTIFYBANKS, BANKFLUX, TRIGGERBANKS,
%   INITIALISEVARIABLES, READMODELINPUTS

% Possible future addition: Should we add some code to check for and
% remove duplicate banks?

Bank.NBanks = sum(Edge.IsBank ~= 0);
Bank.Top = zeros(Bank.NBanks,1);
Bank.Bottom = zeros(Bank.NBanks,1);


%% loop through all banks
BankNo = 0;
EdgeNos = 1:Edge.NEdges;
for ii = EdgeNos(Edge.IsBank~=0) 
    BankNo = BankNo + 1;
    
    DistToBank = abs(Cell.N-Edge.N(ii));
    CellsToBank = [-ii+1:-1,1:(Cell.NCells+1-ii)]';
    
    %% Locate bank top
    if Options.Top == 0
        % Bank top in cell adjacent to IDed bank edge
        Bank.Top(BankNo) = ii - 0.5 - Edge.IsBank(ii)/2;
    else
        % Possible bank top cell if: 
        %  - on the correct side of the bank ID location, and
        %  - within TopCellLim cells of the bank ID location, and
        %  - within TopDistLim distance of bank ID location;
        %  - Or is adjacent cell.
        % (outputs vector of true/false corresponding to each cell)
        PossibleTop = (-Edge.IsBank(ii) * CellsToBank > 0 & ...
                       -Edge.IsBank(ii) * CellsToBank < Options.TopCellLim & ...
                       DistToBank < Options.TopDistLim) | ...
                      (-Edge.IsBank(ii) * CellsToBank == 1);
        switch Options.Top
            case 1
                % NOTE: there are currently no available options for
                % applying a bank top stencil so this is just a holder for
                % potential future implementation.
        end
    end
    
    %% Locate bank toe
    if Options.Bottom == 0 
        % Bank bottom in cell adjacent to IDed bank edge
        Bank.Bottom(BankNo) = ii - 0.5 + Edge.IsBank(ii)/2;
    else
        % Possible toe cell if: 
        %  - on the correct side of the bank ID location, and
        %  - within BotCellLim cells of the bank ID location, and
        %  - within BotDistLim distance of bank ID location;
        %  - Or is adjacent cell.
        % (outputs vector of true/false corresponding to each cell)
        PossibleBottom = (Edge.IsBank(ii) * CellsToBank > 0 & ...
                          Edge.IsBank(ii) * CellsToBank <= Options.BotCellLim & ...
                          DistToBank <= Options.BotDistLim) | ...
                         (Edge.IsBank(ii) * CellsToBank == 1);
        switch Options.Bottom
            case 1 % Max curvature
                Curvature = (Edge.Slope(2:end)-Edge.Slope(1:end-1)) ./ ...
                            (Edge.N(2:end)-Edge.N(1:end-1));
                Curvature(~PossibleBottom) = -999;

                Bank.Bottom(BankNo) = find(Curvature == max(Curvature), 1);
            case 2 % Lowest cell
                Elevation = 9999 * ones(Cell.NCells,1);
                Elevation(PossibleBottom) = Cell.Z(PossibleBottom);
                Bank.Bottom(BankNo) = find(Elevation == min(Elevation), 1);
            case 3 % Max transport cell
                qsTot_flow = zeros(Cell.NCells,1);
                qsTot_flow(PossibleBottom) = Cell.qsTot_flow(PossibleBottom);
                Bank.Bottom(BankNo) = find(qsTot_flow == max(qsTot_flow), 1);
        end
    end
end

end

