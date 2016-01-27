function RMSE = GetModelError(OptIn,OptVar,Inputs,Quiet)
% Wrapper to run model and return error

%% Turn off plots and diagnostics for quicker simulation
if Quiet
    Inputs.Outputs.DiagInt = 0;
    Inputs.Outputs.PlotInt = 0;
    Inputs.Outputs.VideoOut = 0;
end

%% Adjust parameters being optimised
for ii = 1:length(OptVar)
    switch OptVar{ii}
        case 'QsBeRatio'
            Inputs.Opt.Bank.Flux.QsBeRatio = OptIn(ii);
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