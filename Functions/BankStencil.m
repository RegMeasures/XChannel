function [Bank] = BankStencil(Options, Cell, Edge)
% Apply bank erosion stencil to identify top and bottom of bank
%
% Bank.NBanks = Number of active banks currently in model
% Bank.Top = vector (size [NBanks,1]) listing cellNo for the top of each bank
% Bank.Bottom = vector (size [NBanks,1]) listing cellNo for the bottom of each bank

Bank.NBanks = sum(Edge.IsBank ~= 0);
Bank.Top = zeros(Bank.NBanks,1);
Bank.Bottom = zeros(Bank.NBanks,1);


% loop through all banks
BankNo = 0;
EdgeNos = 1:Edge.NEdges;
for ii = EdgeNos(Edge.IsBank~=0) 
    BankNo = BankNo + 1;
    
    DistToBank = Cell.N-Edge.N(ii);
    CellsToBank = [-ii+1:-1,1:(Cell.NCells+1-ii)]';
    
    if Options.Top == 0
        % Bank top in cell adjacent to IDed bank edge
        Bank.Top(BankNo) = ii-0.5 - Edge.IsBank(ii)/2;
    else
        PossibleTop = -Edge.IsBank(ii)*CellsToBank > 0 & ...
                      -Edge.IsBank(ii)*CellsToBank < Options.TopCellLim & ...
                      DistToBank < Options.TopDistLim;
        switch Options.Top
            case 1
                
        end
    end
    
    if Options.Bottom == 0 
        % Bank bottom in cell adjacent to IDed bank edge
        Bank.Bottom(BankNo) = ii-0.5 + Edge.IsBank(ii)/2;
    else
        PossibleBottom = Edge.IsBank(ii)*CellsToBank > 0 & ...
                         Edge.IsBank(ii)*CellsToBank <= Options.BotCellLim & ...
                         DistToBank <= Options.BotDistLim;
        switch Options.Bottom
            case 1 % Max curvature
                
                Curvature = (Edge.Slope(2:end)-Edge.Slope(1:end-1)) ./ ...
                            (Edge.N(2:end)-Edge.N(1:end-1));
                Curvature(~PossibleBottom) = -999;

                Bank.Bottom(BankNo) = find(Curvature == max(Curvature));
        end
    end
end

% do we need to do some cleaning to check for duplicates??

end

