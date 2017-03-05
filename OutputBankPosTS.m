% Post processing routine to loop through all results in outputs folder and
% output time series of bank position.

addpath Functions

OutputsFolder = 'Outputs\';
WL = 212.5;
Radius = 1;

Folders = struct2table(dir(OutputsFolder));
Folders = Folders.name(Folders.isdir);
Folders(ismember(Folders,{'.','..'})) = [];
Scenario = str2double(regexp(Folders,'\d+','match','once'));
[~,SortIndex] = sort(Scenario);
Folders = Folders(SortIndex);

for FolderNo = 1:size(Folders)
    [BankPosition.ModelTime,BankPosition.(Folders{FolderNo,1})] = ...
        BankPosTS(fullfile(OutputsFolder,Folders{FolderNo}, ...
                           ['snapshots_',Folders{FolderNo}]), WL, Radius);
    BankPosition.(Folders{FolderNo,1}) = ...
        BankPosition.(Folders{FolderNo,1}) - ...
        BankPosition.(Folders{FolderNo,1})(1);
end

BankPosition = struct2table(BankPosition);
writetable(BankPosition,'Outputs\BankPosition.csv')