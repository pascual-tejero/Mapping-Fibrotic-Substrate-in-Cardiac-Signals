function [elecIds,elecPos] = extractElecPointIds(elec)

tmp = vtkConnectivityFilter(vtkGenerateIds(elec));
regionIds = unique(tmp.pointData.RegionId);
pointsPerElec = 14;
elecIds = NaN(numel(regionIds), pointsPerElec);
elecPos = NaN(numel(regionIds), 3);
for i = 1:numel(regionIds)
    elecIds(i,:) = tmp.pointData.Ids(tmp.pointData.RegionId==regionIds(i));
    elecPos(i,:) = mean(elec.points(elecIds(i,:),:),1);
end
[elecPos,sortInd] = sortrows(round(1e4*elecPos)/1e4, [2 1]);
elecIds = elecIds(sortInd,:);

end