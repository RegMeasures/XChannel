function [ParValue,ParType] = GetInputParameter(C,ParName,Default,PathName)
% Extract desired parameter from cell array of parameters
% [ParValue,ParType] = GetInputParameter(C,ParName,Default[optional])

RawPar = cell2mat(C{2}(strcmpi(C{1}, ParName)));

if isempty(RawPar)
    % the parameter is missing
    if isempty(Default)
        fprintf('ERROR: Variable %s not found.\n', ParName)
        error('Variable %s not found', ParName)
    end
    ParValue = Default;
    ParType = 0;
    fprintf('Variable %s not found. Set to default value.\n', ParName)
elseif all(ismember(RawPar, '0123456789+-.eEdD'))
    % it's (probably) a number
    ParValue = str2double(RawPar);
    ParType = 1;
    fprintf('%s = %g (numeric)\n',ParName,ParValue)
else
    if ~isempty(PathName)
        FileName = fullfile(PathName,RawPar);
    else
        FileName = RawPar;
    end
    if exist(FileName, 'file')
        % it's a file
        ParValue = csvread(FileName);
        ParType = 2;
        fprintf('%s read in from "%s" (file)\n',ParName,FileName)
    else
        % it's not a number or file so just leave it as a string
        ParValue = RawPar;
        ParType = 3;
        fprintf('%s = %s (string)\n',ParName,ParValue)
    end
end

