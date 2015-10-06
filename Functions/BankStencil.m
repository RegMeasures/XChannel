function [Bank] = BankStencil(Options, Cell, Edge)
% Apply bank erosion stencil to identify top and bottom of bank
%
% Bank.Pos is an matrix with:
%   number of rows = no of cells in cross-section
%   number of cols = no of identified banks (Bank.NBanks)
% for each bank (col) the location of the top and bottom are marked with:
%   -1 = bank top cell
%   1 = bank bottom cell

% create empty BankPos matrix
Bank.NBanks = sum(Edge.IsBank ~= 0);
%Bank.Pos = zeros(Cell.NCells,Bank.NBanks);
Bank.Top = zeros(Bank.NBanks,1);
Bank.Bottom = zeros(Bank.NBanks,1);

% loop through all banks
BankNo = 0;
EdgeNos = 1:Edge.NEdges;
for ii = EdgeNos(Edge.IsBank~=0) 
    BankNo = BankNo + 1;
    switch Options.Top
        case 0 % Bank top in cell adjacent to IDed bank edge
            %Bank.Pos(ii-0.5 - Edge.IsBank(ii)/2,BankNo) = -1;
            Bank.Top(BankNo) = ii-0.5 - Edge.IsBank(ii)/2;
        case 1
    end

    switch Options.Bottom
        case 0 % Bank bottom in cell adjacent to IDed bank edge
            %Bank.Pos(ii-0.5 + Edge.IsBank(ii)/2,BankNo) = 1;
            Bank.Bottom(BankNo) = ii-0.5 + Edge.IsBank(ii)/2;
        case 1
    end
end

% do we need to do some cleaning to check for duplicates??

end

