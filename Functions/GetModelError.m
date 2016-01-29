function RMSE = GetModelError(OptIn,OptVar,Inputs,Quiet)
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
%
% Parameters enabled for optimisation are currently restricted to:
%  - ThetSD
%  - QsBeRatio
%  - BErodibility

%% Turn off plots and diagnostics for quicker simulation
if Quiet
    Inputs.Outputs.DiagInt = 0;
    Inputs.Outputs.PlotInt = 0;
    Inputs.Outputs.VideoOut = 0;
end

%% Adjust parameters being optimised
for ii = 1:length(OptVar)
    switch OptVar{ii}
        case 'Repose'
            Inputs.Opt.Bank.Flux.Repose = OptIn(ii);
        case 'ThetSD'
            Inputs.Opt.Bank.Flux.ThetSD = OptIn(ii);
        case 'QsBeRatio'
            Inputs.Opt.Bank.Flux.QsBeRatio = OptIn(ii);
        case 'BErodibility'
            Inputs.Opt.Bank.Flux.BErodibility = OptIn(ii);
    end
end

%% Run the model and calculate the error
[FinalXS] = XChannelModel(Inputs);
RMSE = rmse(Inputs.Hyd.InitialGeometry(:,3), FinalXS(:,2));
end

function RMSE = rmse(Observation, Model)
% Computed RMSE from comparison of observed and modelled data

Error = Model - Observation;
RMSE = sqrt(mean(Error.^2));

end

function RBE = RightBankError(Observation, Model)
end