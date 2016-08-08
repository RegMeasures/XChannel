function [AbsError, ErrorSign] = GetModelError(OptIn,OptVar,Inputs,BankTestWL,Quiet,LineH)
% Wrapper to run XChannelModel in optimisation routine and return error
%
% RMSE = GetModelError(OptIn,OptVar,Inputs,Quiet)
%
%   where:
% OptIn  = vector of parameter value corresponding to parameters named in
%          OptVar
% OptVar = Parametr names for optimisation (cell array of strings)
% Inputs = XChannelModel inputs as read by ReadModelInputs
% Quiet  = Boolean - True = no outputs, False = outputs as specified in
%          Inputs.Outputs
% LineH  = LineHandle for series of points to be updated with convergance
%          (optional - if omitted no data will be plotted)
%
% Parameters enabled for optimisation are currently restricted to:
%  - Repose
%  - ThetSD
%  - QsBeRatio
%  - BErodibility

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
[FinalXS, WL] = XChannelModel(Inputs);
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