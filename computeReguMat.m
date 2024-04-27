addpath('/Volumes/bordeaux/IBT/matlab/common/vtkToolbox/MATLAB')
addpath(genpath('functions'));

%% read meshes

tissue = vtkRead('../elecTissDist2.1/geometries/tissue.ply');

L = Laplacian_transmural(tissue);

save('../elecTissDist2.1/matrices/L.mat', 'L', '-v7.3');
