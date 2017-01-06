function [EdgeVals] = Center2Edge(CenterVals, C2E_Weights, ActiveEdgeFlag)
%CENTRE2EDGE   Convert cell center parameter to cell edge parameter
%
%   [EdgeVals] = Center2Edge(CenterVals, C2E_Weights, ActiveEdgeFlag)
%   Convert cell centered parameter (CenterVals) to cell edge parameter
%   (EdgeVals) using cell weightings (C2E_Weights) to assign values from
%   the cell to the left or right of each cell (or a combination of the
%   two).
%   
%   Inputs:
%      CenterVals     = NCells x NParameters matrix of cell centre 
%                       parameter values.
%      C2E_Weights    = NCells-1 x 1 matrix of cell weightings. C2E_Weight
%                       values specify the proportion of each edges value
%                       which is derived from the cell to its left.
%                       C2E_Weights are used to control upwind vs central
%                       approaches to calculating edge values.
%      ActiveEdgeFlag = NCells+1 x 1 boolean matrix (mask) identifying 
%                       active cell edges. Non active edges are assigned a
%                       value of 0.
%      
%   Outputs:
%      EdgeVals       = NCells+1 x NParameters matrix of cell edge
%                       parameter values.
%   
%   See also: XCHANNEL.

NParameters = size(CenterVals,2);
C2E_Weights = C2E_Weights*ones(1, NParameters);
EdgeVals = [zeros(1,NParameters);
            CenterVals(1:end-1,:) .* C2E_Weights + ...
            CenterVals(2:end,:) .* (1-C2E_Weights); 
            zeros(1,NParameters)];
EdgeVals(~ActiveEdgeFlag,:) = 0;

end

