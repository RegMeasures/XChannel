function [ParValue,ParType] = GetInputParameter(C,ParName,varargin)
% Extract desired parameter from cell array of parameters
% [ParValue,ParType] = GetInputParameter(C,ParName,Default[optional])

RawPar = cell2mat(C{2}(strcmpi(C{1}, ParName)));

if isempty(RawPar)
    % the parameter is missing
    if isempty(varargin)
        fprintf('ERROR: Variable %s not found.\n', ParName)
        error('Variable %s not found', ParName)
    end
    ParValue = varargin{1};
    ParType = 0;
    fprintf('Variable %s not found. Set to default value.\n', ParName)
elseif all(ismember(RawPar, '0123456789+-.eEdD'))
    % it's (probably) a number
    ParValue = str2double(RawPar);
    ParType = 1;
    fprintf('%s = %g (numeric)\n',ParName,ParValue)
elseif exist(RawPar, 'file')
    % it's a file
    ParValue = csvread(RawPar);
    ParType = 2;
    fprintf('%s read in from "%s" (file)\n',ParName,RawPar)
else
    % it's not a number or file so just leave it as a string
    ParValue = RawPar;
    ParType = 3;
    fprintf('%s = %s (string)\n',ParName,ParValue)
end

