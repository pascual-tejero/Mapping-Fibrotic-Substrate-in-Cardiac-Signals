%% addpath
addpath('/Volumes/bordeaux/IBT10.10/matlab/common/vtkToolbox/MATLAB')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/inverse_problem/BEM')
addpath('/Volumes/bordeaux/IBT10.10/matlab/thirdparty/gptoolbox')
addpath('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/BEM_PrimaryCurrents/LaplacianBlurDownsampling')

%% Initial points
xo=0:500:50000;
yo=0:500:50000;
zo=100:500:3000;

%% meshgrid
[x, y, z]=meshgrid(xo,yo,zo);

%% alpahShape
shp=alphaShape(x(:),y(:),z(:));
plot(shp);

%% boundaryFacets
[bf,P] = boundaryFacets(shp);

%% triangulation
TR = triangulation(bf,P);
tissue_coarse = triangulationToVtk(TR);

%% Read pts electodos
 FileId=fopen('/Volumes/Daten/Benutzer/pt732/Experiments/02_fibrotic_tissue/tissue_ints_10_10/meshes_startstate/2021-05-09_kKPsoUwejU/block.pts');
 npoints=textscan(FileId,'%f',1,'HeaderLines',1);

 points=textscan(FileId,'%f',npoints{2},'MultipleDelimsAsOne',1,'Headerlines',1);
 % now you have the values you want you can put them in a matrix or any variable
 Y=cell2mat(points);
 
%% Electrodes
electrodos= vtkCreateIcosphere(numdiv,radio,centro)


%% load transfer matrix
computeTransferMat.m %load('mwk03_geo/mwk03_Atmv.mat');

%% load source distribution
tmv_fine_vol_point = igb_data('/Volumes/Daten/Benutzer/pt732/Experiments/02_fibrotic_tissue/tissue_ints_10_10/2021-05-09_basic_2001.0_startstate/vm.igb', 'all', 'all');
%% create volume mesh of downsampled fine geometry
heart_coarse_vol_tmp = tetrahedralizeTriangleMesh(heart_coarse_surface);
%% create NablaVm to compute the primary impressed currents J_p = - sigma*nabla*vm 
