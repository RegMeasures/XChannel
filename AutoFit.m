function x = AutoFit(FileName, OptVar, lb, ub)
% Auto fit XChannelModel parameters to achieve best calibration
%
% x = AutoFit(FileName, OptVar, lb, ub)
%
%    where:
% FileName = name of model input file for XChannelModel.m
%            (a final cross-section geometry must be included to fit to)
% OptVar   = list of parametrs to optimise (cell array of strings)
%            (only a limited subset can be optimised: see GetModelError for
%            list of parameters for which optimisation has been enabled)
% lb       = vector of lower bounds corresponding to each named parameter
%            in OptVar
% ub       = vector of upper bounds corresponding to each named parameter
%            in OptVar
% x        = vector of optimised parameter values corresponding to each
%            named parameter in OptVar

addpath('Functions')

%% Read main input settings from model file

% Get file name and path if not specified as input
if ~exist('FileName','var')
    [FileName,FilePath] = uigetfile('*.txt','Select the model input file');
    if isequal(FileName,0)
        error('User selected Cancel')
    end
    FileName = fullfile(FilePath,FileName);
end

% Read model data and options from intput file
[Inputs] = ReadModelInputs(FileName);

%% Do the optimisation...

% List Initial values
x0 = [Inputs.Opt.Bank.Flux.QsBeRatio];

% define function to optimise
fun = @(x)GetModelError(x,OptVar,Inputs,true);

% set optimisation options
options = optimoptions('fmincon');                % default options
options = optimoptions(options,'Display', 'iter');% diagnostic options
options = optimoptions(options,'MaxIter', 10);    % max number of iterations
options = optimoptions(options,'TolX', 0.0001);   % parameter tolerance

% run the optimisation
[x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

%% Plot the final fit
GetModelError(x,OptVar,Inputs,false);

end

