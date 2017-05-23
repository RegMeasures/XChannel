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
[~, ~, raw] = xlsread('Inputs\Scenarios.xlsx','Scenarios','A3:WU51');
% Allocate imported array to column variable names
clear Scenarios
Scenarios.ID        = raw(:,1);
Scenarios.Run       = cell2mat(raw(:,2));
Scenarios.BankID    = cell2mat(raw(:,3));
Scenarios.BankTop   = cell2mat(raw(:,4));
Scenarios.BankBot   = cell2mat(raw(:,5));
Scenarios.BTrigger  = cell2mat(raw(:,6));
Scenarios.BankFlux  = cell2mat(raw(:,7));
Scenarios.StoredBE  = cell2mat(raw(:,8));
Scenarios.UpwindBedload  = cell2mat(raw(:,9));
Scenarios.SlipRatio = cell2mat(raw(:,10));
Scenarios.ConstBend = cell2mat(raw(:,11));
Scenarios.Radius    = cell2mat(raw(:,12));
Scenarios.Geometry  = raw(:,13);
Scenarios.BankTestWL = cell2mat(raw(:,14));
Scenarios.lb        = cell2mat(raw(:,15));
Scenarios.ub        = cell2mat(raw(:,16));
Scenarios.dT        = cell2mat(raw(:,17));
Scenarios.Vgeometry = raw(:,18);
Scenarios.Vradius   = cell2mat(raw(:,19));
Scenarios.VBankTestWL = cell2mat(raw(:,20));
Scenarios.Comments  = raw(:,21);

Scenarios = struct2table(Scenarios);
Scenarios.ID = strtrim(Scenarios.ID);
Scenarios.Geometry = strtrim(Scenarios.Geometry);
Scenarios.Vgeometry(~cellfun(@(x) any(isnan(x)), Scenarios.Vgeometry)) = ...
    strtrim(Scenarios.Vgeometry(~cellfun(@(x) any(isnan(x)), Scenarios.Vgeometry)));

% Clear temporary variables
clear raw

%% Check if output folder exists and if not create one
if ~exist('Outputs','dir')
    mkdir('Outputs')
end

%% Run each selected scenario to optimise bank erosion coefficient

Scenarios.BankCoef = nan(size(Scenarios,1),1);
Scenarios.FinalError = nan(size(Scenarios,1),1);
Scenarios.ValidationError = nan(size(Scenarios,1),1);
Scenarios.SimulationDate = cell(size(Scenarios,1),1);
for ScenNo = 1:size(Scenarios,1)
    if Scenarios.Run(ScenNo)
        % Create folder for scenario outputs
        if ~exist(sprintf('Outputs\\Scenario%s',Scenarios.ID{ScenNo}),'dir')
            mkdir(sprintf('Outputs\\Scenario%s',Scenarios.ID{ScenNo}))
        end

        
        % Set FileName for output
        Inputs.FileName = sprintf('Outputs\\Scenario%s\\Scenario%s.txt',Scenarios.ID{ScenNo},Scenarios.ID{ScenNo});
        
        % Set options
        Inputs.Bank.ID.Approach      = Scenarios.BankID(ScenNo);
        Inputs.Bank.Stencil.Top      = Scenarios.BankTop(ScenNo);
        Inputs.Bank.Stencil.Bottom   = Scenarios.BankBot(ScenNo);
        Inputs.Bank.Trigger.BTrigger = Scenarios.BTrigger(ScenNo);
        Inputs.Bank.Flux.Approach    = Scenarios.BankFlux(ScenNo);
        Inputs.Bank.Update.StoredBE  = Scenarios.StoredBE(ScenNo);
        Inputs.ST.UpwindBedload      = Scenarios.UpwindBedload(ScenNo);
        Inputs.Bank.Flux.SlipRatio   = Scenarios.SlipRatio(ScenNo);
        Inputs.Hyd.ConstBend         = Scenarios.ConstBend(ScenNo);
        Inputs.Hyd.Radius            = Scenarios.Radius(ScenNo);
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
        [Scenarios.BankCoef(ScenNo), Scenarios.FinalError(ScenNo), ...
            Scenarios.ValidationError(ScenNo)] = ...
            AutoFit(Inputs, OptVar, x0, lb, ub, Scenarios(ScenNo,:));
        
        % Output optimised scenarios to file
        Scenarios.SimulationDate{ScenNo} = datestr(now,'dd/mm/yyyy HH:MM:SS');
        writetable(Scenarios, 'Outputs\OptimisationResults.csv')
        
    end
end


