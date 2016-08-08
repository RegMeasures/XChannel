function [x,CalibError,ValidError] = AutoFit(Inputs, OptVar, x0, lb, ub, Scenario)
% Auto fit XChannelModel parameters to achieve best calibration
%
% [x,fval,Vfval] = AutoFit(FileName, OptVar, lb, ub)
%
% inputs:
% Inputs   = Model inputs for XChannelModel.m
%            (a final cross-section geometry must be included to fit to)
% OptVar   = list of parametrs to optimise (cell array of strings)
%            (only a limited subset can be optimised: see GetModelError for
%            list of parameters for which optimisation has been enabled)
% lb       = vector of lower bounds corresponding to each named parameter
%            in OptVar
% ub       = vector of upper bounds corresponding to each named parameter
%            in OptVar
% Vradius  = radius for validation run with calibrated parameter (optional)
% Vgeometry= geometry file for validation run (optional - required if
%            Vradius specified
%
% outputs:
% x        = vector of optimised parameter values corresponding to each
%            named parameter in OptVar
% fval     = Error of final fit
% Vfval    = Error of validation run fit

addpath('Functions')

%% Read model data and options from model input file if required
if ~isstruct(Inputs)
    [Inputs] = ReadModelInputs(Inputs);
end

%% Do the calibration

% basic info
[~, ScenarioName, ~] = fileparts(Inputs.FileName);

% define function to optimise
fun = @(x)GetModelError(x,OptVar,Inputs,Scenario.BankTestWL,true);

% set optimisation options
%options = optimoptions('fmincon');                % default options
%options = optimoptions(options,'Display', 'iter');% diagnostic options
%options = optimoptions(options,'MaxIter', 20);    % max number of iterations
%options = optimoptions(options,'TolX', 0.0001);   % parameter tolerance

options = optimset('fminbnd');  
options = optimset(options,'Display', 'iter');    % diagnostic options
options = optimset(options,'MaxIter', 20);        % max number of iterations
options = optimset(options,'TolX', (ub-lb)/1000); % parameter tolerance
options = optimset(options,'PlotFcns',@PlotFit);  % plotting

% run the optimisation
%[x,CalibError] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
[x,CalibError] = fminbnd(fun,lb,ub,options);

%% Plot the final fit (and save plot + animation)
[CalibError,ErrorSign] = GetModelError(x,OptVar,Inputs,Scenario.BankTestWL,false);
CalibError = CalibError * ErrorSign;

%% Validate model
if ~isnan(Scenario.Vradius)
    % model setup
    ValidationInputs = Inputs;
    ValidationInputs.Hyd.Radius = Scenario.Vradius;
    ValidationInputs.Outputs.CsvInt = 99999999; % only output final XS shape for validation run
    if ~isnan(Scenario.Vgeometry{1})
        ValidationInputs.Hyd.InitialGeometry   = csvread(Scenario.Vgeometry{1});
    else
        error('Geometry file %s not found',Scenario.Vgeometry{1})
    end
    ValidationInputs.FileName = [Inputs.FileName(1:end-4),'_Validation.txt'];
    for ii = 1:length(OptVar)
        switch OptVar{ii}
            case 'Repose'
                ValidationInputs.Opt.Bank.Flux.Repose = x(ii);
            case 'ThetSD'
                ValidationInputs.Opt.Bank.Flux.ThetSD = x(ii);
            case 'QsBeRatio'
                ValidationInputs.Opt.Bank.Flux.QsBeRatio = x(ii);
            case 'BErodibility'
                ValidationInputs.Opt.Bank.Flux.BErodibility = x(ii);
        end
    end
    % run model
    [ValidError,ErrorSign] = GetModelError(x,OptVar,ValidationInputs,Scenario.VBankTestWL,false);
    ValidError = ValidError * ErrorSign;
else
    ValidError = nan;
end


%% Nested function for plotting optimisation progress
    function stop = PlotFit(x, optimValues, state)
        % plotfcn for autofit optimisation routine
        stop = false;
        hold on;
        switch state
            case 'init'
                %ylabel('RMSE [m]')
                ylabel('Error in right bank position (m)')
                xlabel(OptVar{1})
                title(ScenarioName)
                xlim([lb, ub])
                ax = gca;
                %ax.YScale = 'log';
            case 'iter'
                plot(x(1),optimValues.fval,'bx');
                drawnow
            case 'interrupt'
                stop = true;
            case 'done'
                plot(x(1),optimValues.fval,'ro');
                drawnow
                saveas(gcf, [Inputs.FileName(1:end-4), '_OptimisationPlot'],'png')
                hold off
        end
    end
end


