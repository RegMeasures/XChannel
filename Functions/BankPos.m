function BankY = BankPos(Dist, Bed, Radius, WL)
%BANKPOS   compute outside bank waters edge position
if Radius > 0
    % Left hand bend, right bank is outside bank
    BankCell = find(Bed<WL, 1, 'last');
    if BankCell < length(Dist)
        BankY = Dist(BankCell) + ...
                      (WL - Bed(BankCell)) * ...
                      (Dist(BankCell+1) - Dist(BankCell)) / ...
                      (Bed(BankCell+1) - Bed(BankCell));
    else
        BankY = Dist(end);
    end
else
    % Right hand bend, left bank is outside bank
    BankCell = find(Bed<WL, 1, 'first');
    if BankCell > 1
        BankY = Dist(BankCell) - ...
                      (WL - Bed(BankCell)) * ...
                      (Dist(BankCell) - Dist(BankCell-1)) / ...
                      (Bed(BankCell-1) - Bed(BankCell));
    else
        BankY = Dist(1);
    end
end    

end