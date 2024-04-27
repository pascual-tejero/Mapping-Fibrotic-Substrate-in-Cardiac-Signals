clc
clear
close all

%% addpath
addpath('/Volumes/bordeaux/IBT10.10/matlab/common/vtkToolbox/MATLAB')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/inverse_problem/BEM')
addpath('/Volumes/bordeaux/IBT10.10/matlab/thirdparty/gptoolbox')
% addpath('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/BEM_PrimaryCurrents/LaplacianBlurDownsampling')
addpath('/Volumes/beaujolais/Benutzer/js191/Pascual/tikhonov')
addpath('/Volumes/beaujolais/Benutzer/js191/Pascual/operators')

%% Initial points
xo=-25000:1000:25000;
yo=-25000:1000:25000;
zo=-500:1000:500;

%% meshgrid
[x, y, z] = meshgrid(xo,yo,zo); %Return a 3D grid, according to the inputs

%% alpahShape
shp = alphaShape(x(:),y(:),z(:)); %C2reates a bounding area or volume that envelops a set of 3-D points
% plot(shp);

%% boundaryFacets
%Boundary facets bf, returned as a matrix. bf is of size m-by-n, where m is 
%the number of boundary facets and n is the number of vertices per facet.
%Vertex coordinates P, returned as a matrix. P is of size N-by-dim, where N 
%is the number of points on the boundary of the alpha shape and dim is either
%2 or 3 (for either a 2-D or 3-D alpha shape)

[bf,P] = boundaryFacets(shp); %Returns a matrix representing the facets that 
%make up the boundary of the alpha shape


%% triangulation
%Use triangulation to create an in-memory representation of any 2-D or 3-D 
%triangulation data that is in matrix format
TR = triangulation(bf,P);

vtk_coarse = triangulationToVtk(TR);

%% Bounding box
xob=-30000:5000:30000;
yob=-30000:5000:30000;
zob=-1000:1000:2000;

[xb, yb, zb] = meshgrid(xob,yob,zob); %Return a 3D grid, according to the inputs

shpb = alphaShape(xb(:),yb(:),zb(:)); %Creates a bounding area or volume that envelops a set of 3-D points

[bfb,Pb] = boundaryFacets(shpb); %Returns a matrix representing the facets that 

TRb = triangulation(bfb,Pb);
vtk_coarse_bbox = triangulationToVtk(TRb);


%% Read pts electodos (unknown)

% xo=-10500:3000:10500;
% yo=-10500:3000:10500;

% xo=-16500:3000:16500;
% yo=-16500:3000:16500;

xo=-22500:3000:22500;
yo=-22500:3000:22500;
 %% meshgrid
[x, y] = meshgrid(xo,yo); %Return a 3D grid, according to the inputs

points=[x(:),y(:),1000.*ones(length(xo)*length(yo),1)]; %distance = 1000 or 1500 
 
%% Electrodes
%We do not want to have so much subdivisions in the spheres, becuase it will 
%increase the computational cost

%Electrode mesh

[x_sp,y_sp,z_sp] = sphere(4);

r = 1; %Radio

%Creation of the spheres
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
%Boundary facets bf, returned as a matrix. bf is of size m-by-n, where m is 
%the number of boundary facets and n is the number of vertices per facet.
%Vertex coordinates P, returned as a matrix. P is of size N-by-dim, where N 
%is the number of points on the boundary of the alpha shape and dim is either
%2 or 3 (for either a 2-D or 3-D alpha shape).

[bf,P] = boundaryFacets(shp_elec); %Returns a matrix representing the facets that 
%make up the boundary of the alpha shape

%% triangulation
%Use triangulation to create an in-memory representation of any 2-D or 3-D 
%triangulation data that is in matrix format
TR_elec = triangulation(bf,P);

vtk_electrodes = triangulationToVtk(TR_elec); %vtkfile

vtkWrite(vtk_coarse,'sources.vtk',false,'ascii');
vtkWrite(vtk_electrodes,'electrodes.vtk',false,'ascii');

%% Transfer matrix
% exist_A = exist('A');
% 
% if exist_A == 0
    A = computeTransferMat_A2(vtk_coarse, vtk_electrodes,vtk_coarse_bbox);  %This needs to be calculated before hand TODO create a function to generate the A matrix if not exists
% else
%     load('A.mat');
% end
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/A.mat', 'A', '-v7.3');

%% Matrix L (regularization term)
L = Laplacian_transmural(vtk_coarse);

% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/L.mat', 'L', '-v7.3');
%% Matrix B

%Load electrogram data
phie = igb_data('/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02_6/d1000_e16x16/2021-08-19_basic_3000.0/phie_recovery.igb', 'all', 'all');

B = phie;
% B = B(:,end-300:end);

%Concatenation of a matrix of zeros and the matrix phie
% B = cat(2,zeros(size(phie,1),5),phie);

% refSig = zeros(size(B,1),1); %mean(phie_t(:,:),1)';

% B = B - repmat(mean(B,1),size(B,1),1);
% B = B - repmat(mean(B,2),1,size(B,2));

save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/B.mat', 'B', '-v7.3');


%% SVD
[U,Q,S,M] = tikhonov_gsvd_decompose(A,L);
lambda = tikhonov_gsvd_Lcurve(B,U,Q,S,M,[-5 2]);
[X,~,~] = tikhonov_gsvd(B,U,Q,S,M,lambda);
% [X,resNorm,solNorm] = tikhonov(B, A, L, 1e-5);

%%
Xm = movmean(X,[3 3],2);

%%
sigma = 1;
sig = gaussFiltfilt(Xm, sigma);
derivSig = gradient(sig); %derivSig = (sig(:,3:end)-sig(:,1:end-2))/2; %
[~,at] = min(derivSig,[],2);

comb = vtkAppendPolyData({vtk_coarse, vtk_electrodes});

qtrip(comb.points, comb.cells, Xm, 1041, [0 90 90], 5, [min(min(Xm)) max(max(Xm))], 10);
pause(0.5)
qtrip(comb.points, comb.cells, at, 1042, [0 90 90], 5, [0 max(at)], 1);

save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/X.mat', 'X', '-v7.3');
