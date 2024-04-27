clc
clear
close all

%% addpath
addpath('/Volumes/bordeaux/IBT10.10/matlab/common/vtkToolbox/MATLAB')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/inverse_problem/BEM')
addpath('/Volumes/bordeaux/IBT10.10/matlab/thirdparty/gptoolbox')
%addpath('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/BEM_PrimaryCurrents/LaplacianBlurDownsampling')
addpath('/Volumes/beaujolais/Benutzer/js191/Pascual/tikhonov')
addpath('/Volumes/beaujolais/Benutzer/js191/Pascual/operators')
%% Initial points
xo=-25000:1000:25000;
yo=-25000:1000:25000;
zo=-500:1000:500;

%% meshgrid
[x, y, z] = meshgrid(xo,yo,zo); %Return a 3D grid, according to the inputs

%% alpahShape
shp = alphaShape(x(:),y(:),z(:)); %Creates a bounding area or volume that envelops a set of 3-D points

%% boundaryFacets
[bf,P] = boundaryFacets(shp); %Returns a matrix representing the facets that 
%make up the boundary of the alpha shape

%Boundary facets bf, returned as a matrix. bf is of size m-by-n, where m is 
%the number of boundary facets and n is the number of vertices per facet.
%Vertex coordinates P, returned as a matrix. P is of size N-by-dim, where N 
%is the number of points on the boundary of the alpha shape and dim is either
%2 or 3 (for either a 2-D or 3-D alpha shape).
%% triangulation
%Use triangulation to create an in-memory representation of any 2-D or 3-D 
%triangulation data that is in matrix format
TR = triangulation(bf,P);
vtk_coarse = triangulationToVtk(TR);

%% Read pts electodos (unknown)

xo=-10500:3000:10500;%-3500:1000:3500;
yo=-10500:3000:10500;%-3500:1000:3500;

 %% meshgrid
[x, y] = meshgrid(xo,yo); %Return a 3D grid, according to the inputs

points=[x(:),y(:),500.*ones(64,1)];
 
%% Electrodes
%We do not want to have so much subdivisions, becuase it will increase the
%computational cost

%electrodes mesh

[x_sp,y_sp,z_sp] = sphere(4);

r = 10; %Radio

Xr = x_sp * r;
Yr = y_sp * r;
Zr = z_sp * r;

Xr = Xr(:);
Yr = Yr(:);
Zr = Zr(:);

electrode_sphere = [];

for i=1:size(points,1)
    X_sphere = points(i,1) + Xr;
    Y_sphere = points(i,2) + Yr;
    Z_sphere = points(i,3) + Zr;

    temp = [X_sphere Y_sphere Z_sphere];
    electrode_sphere = [electrode_sphere; temp];

end

shp_elec = alphaShape(electrode_sphere(:,1),electrode_sphere(:,2),electrode_sphere(:,3));
% plot(shp_elec);

%% boundaryFacets
[bf,P] = boundaryFacets(shp_elec); %Returns a matrix representing the facets that 
%make up the boundary of the alpha shape

%Boundary facets bf, returned as a matrix. bf is of size m-by-n, where m is 
%the number of boundary facets and n is the number of vertices per facet.
%Vertex coordinates P, returned as a matrix. P is of size N-by-dim, where N 
%is the number of points on the boundary of the alpha shape and dim is either
%2 or 3 (for either a 2-D or 3-D alpha shape).
%% triangulation
%Use triangulation to create an in-memory representation of any 2-D or 3-D 
%triangulation data that is in matrix format
TR_elec = triangulation(bf,P);

vtk_electrodes = triangulationToVtk(TR_elec);

vtkWrite(vtk_coarse,'sources.vtk',false,'ascii');
vtkWrite(vtk_electrodes,'electrodes.vtk',false,'ascii');

%% Transfer matrix

A = computeTransferMat_A(vtk_coarse, vtk_electrodes);  %This needs to be calculated before hand TODO create a function to generate the A matrix if not exists


%% Matrix L (regularization term)
L = Laplacian_transmural(vtk_coarse);

%% Matrix B

%B = igb_data('/Volumes/beaujolais/Benutzer/pt732/phie_recovery.igb', 'all', 'all');%/Volumes/Daten/pt732/2021-05-28_basic_1000.0/phie_recovery.igb', 'all', 'all');

%B = igb_data('/Volumes/beaujolais/Benutzer/pt732/phie_recovery_r1.igb', 'all', 'all');

B = igb_data('/Volumes/beaujolais/Benutzer/pt732/phie_recovery_r1_d1000_e8x8.igb', 'all', 'all');


%B = cat(2,zeros(size(phie,1),5),phie);

refSig = zeros(size(B,1),1);

B = B - repmat(refSig, 1, size(B,2));

%% Reconstruct
%lambda = 1e-3;
%%lambda = tikhonovLcurve_corner(B, A, L, true, [-6 3], 50, 1, gca);
%X = tikhonov_direct(B, A, L, lambda, false);

%% SVD
[U,Q,S,M] = tikhonov_gsvd_decompose(A,L);
lambda = tikhonov_gsvd_Lcurve(B,U,Q,S,M,[-5 2]);
[X,~,~] = tikhonov_gsvd(B,U,Q,S,M,lambda);

%%
Xm = movmean(X,[3 3],2);

%%
sigma = 1;
sig = gaussFiltfilt(Xm, sigma);
derivSig = gradient(sig);
[~,at] = min(derivSig,[],2);

comb = vtkAppendPolyData({vtk_coarse, vtk_electrodes});

qtrip(comb.points, comb.cells, Xm, 1041, [0 90 90], 0.5, [min(min(Xm)) max(max(Xm))], 10);
pause(0.5)
qtrip(comb.points, comb.cells, at, 1042, [0 90 90], 1, [0 max(at)], 1);
