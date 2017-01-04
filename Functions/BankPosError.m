function [AbsBankError, ErrorSign] = BankPosError(Dist, ObsBed, ...
                                                  ModelBed, Radius, WL)
%BANKPOSERROR   Error (diff) in bank position between two cross-sections
%Calculate the horizontal difference in outside bend bank position between
%two cross-sections. Bank position is defined as the point at which a
%horizontal water surface (defined by WL) intersects the outside bank. Bank
%position is linearly interpolated between cell-center positions.
%
%   [AbsBankError, ErrorSign] = BANKPOSERROR(Dist, ObsBed, ModelBed, ...
%                                            Radius, WL)
%   Inputs:
%      Dist     = vector of cross-channel positions of cell-centers 
%                 (in increasing order)
%      ObsBed   = surveyed bed level of cell centers 
%                 (same size array as Dist)
%      ModelBed = modelled bed level of cell centers
%      Radius   = bend radius (used only to determin which side of the XS 
%                 is the outside of the bend)
%      WL = Water level (used for computing bank position)
%   
%   Outputs:
%      AbsBankError = Absoulte error in bank position calculated as difference
%                     between waters edge position for observed and modelled
%                     final topographies
%      ErrorSign = Direction of error: +ve = too much bank erosion
%                                      -ve = not enough bank erosion

% Validate inputs
assert(isvector(Dist))
assert(isvector(ObsBed))
assert(isvector(ModelBed))
assert(isequal(numel(Dist),numel(ObsBed),numel(ModelBed)))

% Calculate the error in bank position
ModelBank = BankPos(Dist, ModelBed, Radius, WL);
ObsBank = BankPos(Dist, ObsBed, Radius, WL);
BankError = ModelBank - ObsBank;
if Radius > 0
    % Left hand bend, right bank is outside bank
    ErrorSign = sign(BankError);
else
    % Right hand bend, left bank is outside bank
    ErrorSign = -sign(BankError);
end
AbsBankError = abs(BankError);
end

function BankY = BankPos(Dist, Bed, Radius, WL)
% Private function to compute outside bank waters edge position
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