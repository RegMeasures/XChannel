function [ParValue,ParType] = GetInputParameter(C, ParName, Default, ...
                                                PathName)
%GETINPUTPARAMETER   Parse cell array to extract named input parameter
%GETINPUTPARAMETER parses a cell array and extracts the value associated
%with a specific named parameter. Can handle parameters which are numeric, 
%strings, missing, or are filenames of csv formatted numeric files.
%
%   [ParValue, ParType] = GetInputParameter(C, ParName, ...
%                                           Default(optional), ...
%                                           PathName(optional))
%
%   Inputs:
%      C        = Cell array of parameter names and values where:
%                    C(1,1) = [Nx1] cell array of parameter names
%                    C(1,2) = [Nx1] cell array of parameter values
%                 Typically C will have been read in from a text file using 
%                 textscan.
%      ParName  = Name of desired parameter to extract (string). Note that
%                 parameter names are not case sensitive.
%      Default  = Default parameter value if not found (optional).
%      Pathname = Path name to location of any filename variables
%                 (optional).
%
%   Outputs depend on parameter type:
%      Missing parameters
%         If a default is provided then: ParValue = Default; ParType = 0;
%         If no default then GETINPUTPARAMETER returns an error
%      Numeric parameters
%         If the parameter value can be interpretted as a number then:
%         ParValue = double; ParType = 1;
%      Csv file parameters
%         If a file exists with name = parameter value then it is assumed
%         the parameter refers to a csv file.
%         ParValue = matrix of numeric values read in from CSV file using
%         csvread; ParType = 2;
%      String parameter
%         All parameters which exist, are non-numeric and do not refer to 
%         a existing file are read in as strings: 
%         ParValue = string; ParType = 3;
%
%   See also: READMODELINPUTS TEXTSCAN CSVREAD

RawPar = cell2mat(C{1,2}(strcmpi(C{1,1}, ParName)));

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

