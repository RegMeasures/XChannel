function [Times,BankCoord] = BankPosTS(SnapshotFolder, WL, Radius)
%BANKPOSTS   Extract timeseries of bank position from XChannel snapshots

% get list of snapshot files
Snapshots = dir(fullfile(SnapshotFolder,'*.out'));
Snapshots = struct2table(Snapshots);
Snapshots = Snapshots.name;
NoOfFiles = size(Snapshots,1);

% extract times from filename
Times = regexp(Snapshots,'\d+','match','once');
Times = str2double(Times);

% Loop through files calculating bank position
BankCoord = nan(NoOfFiles,1);
for FileNo = 1:NoOfFiles
    Geometry = csvread(fullfile(SnapshotFolder,Snapshots{FileNo}));
    BankCoord(FileNo) = BankPos(Geometry(:,1), Geometry(:,2), Radius, WL);
end

[Times,SortIndex] = sort(Times);
BankCoord = BankCoord(SortIndex);
end

