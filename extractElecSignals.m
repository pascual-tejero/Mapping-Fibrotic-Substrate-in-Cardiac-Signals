addpath('/Volumes/bordeaux/IBT/matlab/common/vtkToolbox/MATLAB')

%%
vol_all = vtkRead('/Volumes/beaujolais/Benutzer/js191/Steffen/mesh/mesh_scn3_elecTissDist2.1.vtk');

% translate mesh so that tissue center is at [0,0,0] and scale to mm
vol_tiss = vtkThreshold(vol_all, 'cells', 'elemTag', [1 1]);
vol_all.points = 1e-3 * (vol_all.points-(max(vol_tiss.points,[],1)+min(vol_tiss.points,[],1))/2);

elec = vtkRead('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/intracardiac_inverse/elecTissDist2.1/geometries/electrodes.ply');
elec = vtkMapIds(vol_all, elec);
elecIdsInVol = elec.pointData.MappedIds;

elecIdsInSur = extractElecPointIds(elec);

elecIds = elecIdsInVol(elecIdsInSur);
refIds = find(abs(vol_all.points(:,3)-max(vol_all.points(:,3))) < 1e-4);

load('/Volumes/beaujolais/Benutzer/js191/Steffen/simulation_2019-12-24_scn3_elecTissDist_2.1mm/phie.mat');

B = reshape(phie(:,elecIds+1), size(phie,1), size(elecIds,1), size(elecIds,2));
B = mean(B,3)';

refSig = mean(phie(:,refIds+1),2)';
B = B - repmat(refSig, size(B,1), 1);

% save('../elecTissDist2.1/signals/B.mat', 'B', '-v7.3');
