addpath('/Volumes/bordeaux/IBT/matlab/common/vtkToolbox/MATLAB')
addpath(genpath('functions'));

%% read meshes

tiss = vtkRead('../elecTissDist2.1/geometries/tissue.ply');
elec = vtkRead('../elecTissDist2.1/geometries/electrodes.ply');
bath = vtkRead('../elecTissDist2.1/geometries/bath.ply');

tiss_mesh = PrepareTriangleMesh(tiss.points, tiss.cells);
elec_mesh = PrepareTriangleMesh(elec.points, elec.cells);
bath_mesh = PrepareTriangleMesh(bath.points, bath.cells);

%% extract point ids of electrodes

elecIds = extractElecPointIds(elec);

%% extract point ids of top of bath for defining the reference level

refIds = find(abs(bath.points(:,3)-max(bath.points(:,3))) < 1e-4)';

%% compute transfer matrix

c_in = 0.05; % intracellular conductivity of source domain

% bulk conductivities
ci = [0.2, 1e3, 0.7]; % inside of compartments
co = [0.7, 0.7, 0.0]; % outside of compartments

meshes = {tiss_mesh, elec_mesh, bath_mesh}; % source surface has to come first
electrodesSurface = 2;
refSurface = 3;

zerolevel = [refSurface refIds];
[Ltmv,inds] = TransferMatrix_TMV_Linear(meshes, ci, co, c_in, zerolevel);

A_allPoints = Ltmv(inds{electrodesSurface},:);

A = NaN(size(elecIds,1), size(A_allPoints,2));
for i = 1:size(elecIds,1)
    A(i,:) = mean(A_allPoints(elecIds(i,:),:), 1);
end

% save('../elecTissDist2.1/matrices/A.mat', 'A', '-v7.3');
