function [AbsError, ErrorSign] = GetModelError(OptIn, OptVar, Inputs, ...
                                               BankTestWL, Quiet, LineH)
%GETMODELERROR   Run XChannel with modified parameters & return error
%
%   [AbsError, ErrorSign] = GetModelError(OptIn, OptVar, Inputs, ...
%                                         BankTestWL ,Quiet, LineH)
%
%   Inputs:
%      OptIn  = vector of parameter value corresponding to parameters named
%               in OptVar
%      OptVar = Parameter names for optimisation (cell array of strings). 
%               Parameters enabled for optimisation are restricted to:
%                  - Repose
%                  - ThetSD
%                  - QsBeRatio
%                  - BErodibility
%      Inputs = XChannel inputs as read by ReadModelInputs (struct)
%      Quiet  = True = no outputs, False = outputs as specified by 
%               Inputs.Outputs (boolean)
%      LineH  = LineHandle for series of points to be updated with
%               convergance (optional - if omitted no data will be plotted)
%   
%   Outputs:
%      AbsError  = Absoulte error in horizontal bank position calculated as 
%                  difference between waters edge position for observed 
%                  and modelled final cross-section
%      ErrorSign = Direction of error: +ve = too much bank erosion
%                                      -ve = not enough bank erosion
%
%   See also: XCHANNEL, BANKPOSERROR, AUTOFIT, READMODELINPUTS.

%% Turn off plots and diagnostics for quicker simulation
if Quiet
    Inputs.Outputs.DiagInt = 0;
    Inputs.Outputs.PlotInt = 0;
    Inputs.Outputs.VideoOut = 0;
    Inputs.Outputs.CsvInt = 0;
end

%% Adjust parameters being optimised
for ii = 1:length(OptVar)
    switch OptVar{ii}
        case 'Repose'
            Inputs.Bank.Flux.Repose = OptIn(ii);
        case 'ThetSD'
            Inputs.Bank.Flux.ThetSD = OptIn(ii);
        case 'QsBeRatio'
            Inputs.Bank.Flux.QsBeRatio = OptIn(ii);
        case 'BErodibility'
            Inputs.Bank.Flux.BErodibility = OptIn(ii);
    end
end

%% Run the model and calculate the error
[FinalXS, ~] = XChannel(Inputs);
%AbsError = rmse(Inputs.Hyd.InitialGeometry(:,3), FinalXS(:,2));
[AbsError, ErrorSign] = BankPosError(FinalXS(:,1), ...
                                     Inputs.Hyd.InitialGeometry(:,3), ...
                                     FinalXS(:,2), Inputs.Hyd.Radius, ...
                                     BankTestWL);

%% Update the optimisation figure (optional)
if exist('LineH','var')
    set(LineH,'XData',[LineH.XData, OptIn(1)],...
              'YData',[LineH.YData, AbsError*ErrorSign])
    drawnow
end
end

function RMSE = rmse(Observation, Model)
% Computed RMSE from comparison of observed and modelled data

BedError = Model - Observation;
RMSE = sqrt(mean(BedError.^2));
end