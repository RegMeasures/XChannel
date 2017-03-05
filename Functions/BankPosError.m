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
%      AbsBankError = Absoulte error in bank position calculated as 
%                     difference between waters edge position for observed 
%                     and modelled final cross-section
%      ErrorSign = Direction of error: +ve = too much bank erosion
%                                      -ve = not enough bank erosion
%
%   See also: GETMODELERROR, XCHANNEL.

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