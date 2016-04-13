function [x,fval,Vfval] = AutoFit(Inputs, OptVar, x0, lb, ub, Vradius, Vgeometry)
% Auto fit XChannelModel parameters to achieve best calibration
%
% [x,fval,Vfval] = AutoFit(FileName, OptVar, lb, ub)
%
% inputs:
% FileName = name of model input file for XChannelModel.m
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

%% Do the optimisation...

% basic info
[~, ScenarioName, ~] = fileparts(Inputs.FileName);

% define function to optimise
fun = @(x)GetModelError(x,OptVar,Inputs,true);

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
%[x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
[x,fval] = fminbnd(fun,lb,ub,options);

%% Plot the final fit (and save plot + animation)
GetModelError(x,OptVar,Inputs,false);

%% Validate model
if exist('Vradius','var')
    % model setup
    ValidationInputs = Inputs;
    ValidationInputs.Hyd.Radius = Vradius;
    if exist(Vgeometry, 'file')
        ValidationInputs.Hyd.InitialGeometry   = csvread(Vgeometry);
    else
        error('Geometry file %s not found',Vgeometry)
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
    Vfval = GetModelError(x,OptVar,ValidationInputs,false);
else
    Vfval = nan;
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


