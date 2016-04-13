%% Sensitivity analysis of bank erosion routines
% Loop through each scenario
% Automatically calibrate bank erosion coefficient
% Output table of calibrated coefficients and final fit, animation, still 
% plot of final calibration, and optimisation convergence plot.
% Optimisation based on right bank waters edge position.

addpath('Functions')

%% Load base model data
FileName = 'Inputs\SelwynModel.txt';
[Inputs] = ReadModelInputs(FileName);

%% Load sensitivity scenarios table
% Read excel file
[~, ~, raw] = xlsread('Inputs\Scenarios.xlsx','Scenarios','A3:P66');
% Allocate imported array to column variable names
clear Scenarios
Scenarios.ID        = (1:size(raw,1))';
Scenarios.Run       = cell2mat(raw(:,2));
Scenarios.BankID    = cell2mat(raw(:,3));
Scenarios.BankTop   = cell2mat(raw(:,4));
Scenarios.BankBot   = cell2mat(raw(:,5));
Scenarios.BTrigger  = cell2mat(raw(:,6));
Scenarios.BankFlux  = cell2mat(raw(:,7));
Scenarios.StoredBE  = cell2mat(raw(:,8));
Scenarios.Geometry  = raw(:,9);
Scenarios.Flow      = raw(:,10);
Scenarios.lb        = cell2mat(raw(:,11));
Scenarios.ub        = cell2mat(raw(:,12));
Scenarios.dT        = cell2mat(raw(:,13));
Scenarios.Vgeometry = raw(:,14);
Scenarios.Vradius   = cell2mat(raw(:,15));
Scenarios.Comments = raw(:,16);
Scenarios = struct2table(Scenarios);
% Clear temporary variables
clear raw

%% Run each selected scenario to optimise bank erosion coefficient

Scenarios.BankCoef = nan(size(Scenarios,1),1);
Scenarios.FinalError = nan(size(Scenarios,1),1);
Scenarios.ValidationError = nan(size(Scenarios,1),1);
for ScenNo = 1:size(Scenarios,1)
    if Scenarios.Run(ScenNo)
        % Set FileName for output
        Inputs.FileName = sprintf('Outputs\\Scenario%i.txt',ScenNo);
        
        % Set options
        Inputs.Bank.ID.Approach      = Scenarios.BankID(ScenNo);
        Inputs.Bank.Stencil.Top      = Scenarios.BankTop(ScenNo);
        Inputs.Bank.Stencil.Bottom   = Scenarios.BankBot(ScenNo);
        Inputs.Bank.Trigger.BTrigger = Scenarios.BTrigger(ScenNo);
        Inputs.Bank.Flux.Approach    = Scenarios.BankFlux(ScenNo);
        Inputs.Bank.Update.StoredBE  = Scenarios.StoredBE(ScenNo);
        if exist(Scenarios.Geometry{ScenNo}, 'file')
            Inputs.Hyd.InitialGeometry   = csvread(Scenarios.Geometry{ScenNo});
        else
            error('Geometry file %s not found',Scenarios.Geometry{ScenNo})
        end
        Inputs.Time.dT                   = Scenarios.dT(ScenNo);

        % set parameters to optimise and appropriate range
        switch Inputs.Bank.Flux.Approach
            case 1
                OptVar = {'Repose'};
                x0 = Inputs.Bank.Flux.Repose;
            case 2
                OptVar = {'ThetSD'};
                x0 = Inputs.Bank.Flux.ThetSD;
            case 3
                OptVar = {'QsBeRatio'};
                x0 = Inputs.Bank.Flux.QsBeRatio;
            case 4
                OptVar = {'BErodibility'};
                x0 = Inputs.Bank.Flux.BErodibility;
        end
        lb = Scenarios.lb(ScenNo);
        ub = Scenarios.ub(ScenNo);
        
        % do the optimisation
        if ~isnan(Scenarios.Vradius(ScenNo))
            [Scenarios.BankCoef(ScenNo), Scenarios.FinalError(ScenNo), Scenarios.ValidationError(ScenNo)] = ...
                AutoFit(Inputs, OptVar, x0, lb, ub, ...
                        Scenarios.Vradius(ScenNo), ...
                        Scenarios.Vgeometry{ScenNo});
        else
            [Scenarios.BankCoef(ScenNo), Scenarios.FinalError(ScenNo)] = ...
                AutoFit(Inputs, OptVar, x0, lb, ub);
        end
        
        % Output optimised scenarios to file
        writetable(Scenarios(Scenarios.Run,:), 'Outputs\OptimisationResults.csv')
        
    end
end


