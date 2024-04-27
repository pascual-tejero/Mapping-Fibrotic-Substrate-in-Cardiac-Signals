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
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/EGMDataProcessing/FilteringAndProcessing')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/EGMDataProcessing/activityAndDPs')
addpath('activationtimes-master/computeGradientAndHessian')

%% Initial points
xo_points=-25000:1000:25000; % x axis
yo_points=-25000:1000:25000; % y axis
zo_points=-500:1000:500; % z axis

%% meshgrid
[x, y, z] = meshgrid(xo_points,yo_points,zo_points); %Return a 3D grid, according to the inputs

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

points=[x(:),y(:),1500.*ones(length(xo)*length(yo),1)]; %distance = 1000 or 1500 
 
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

% if exist_A == 0
    A = computeTransferMat_A(vtk_coarse, vtk_electrodes,vtk_coarse_bbox);  %This needs to be calculated before hand TODO create a function to generate the A matrix if not exists
% else
%     load('A.mat');
% end
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/A.mat', 'A', '-v7.3');

%% Matrix L (regularization term)
L = Laplacian_transmural(vtk_coarse);

% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/L.mat', 'L', '-v7.3');
%% Matrix B

%Load electrogram data
phie = igb_data('/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1500_e16x16/2021-06-11_basic_1000.0/phie_recovery.igb', 'all', 'all');

B = phie;
% B = B(:,1:end);
B = B(:,1:125);
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

qtrip(comb.points, comb.cells, Xm, 1041, [0 90 90], 0.5, [min(min(Xm)) max(max(Xm))], 10);
pause(0.5)
qtrip(comb.points, comb.cells, at, 1042, [0 90 90], 5, [0 max(at)], 1);

save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/X.mat', 'X', '-v7.3');


%%
% V = double(vtk_coarse.points); F = double(vtk_coarse.cells);
% G = grad(V,F);
% W = vtk_coarse;
% NablaVm = squeeze(reshape(G*W,size(F,1),size(V,2),size(W,2)));
% 
% NablaVm_coarse = NablaVm;
% 
% % TODO: apply second iteration of laplacian blurring to NablaVm!
% %% calculate distance functions, 
% 
% % make cellData to pointData?
% cellCenters_struct = vtkCellCentroids(heart_coarse_vol_tmp);
% cellCenters_w_fib = cellCenters_struct.points;
% clearvars cellCenters_struct
% D = pdist2(single(cellCenters_w_fib),single(heart_coarse_surface.points));
% D_torso = pdist2(single(cellCenters_w_fib),single(torso_coarse_surface.points));
% 
% %% calculate element volumes
% vol_elem = zeros(size(heart_coarse_vol_tmp.cells,1),1);
% for i = 1:size(heart_coarse_vol_tmp.cells,1)
%     
%     xyz = [1 heart_coarse_vol_tmp.points(heart_coarse_vol_tmp.cells(i,1),:) ; 1 heart_coarse_vol_tmp.points(heart_coarse_vol_tmp.cells(i,2),:) ;...
%         1 heart_coarse_vol_tmp.points(heart_coarse_vol_tmp.cells(i,3),:) ; 1 heart_coarse_vol_tmp.points(heart_coarse_vol_tmp.cells(i,4),:)]; 
%     vol_elem(i,1) = det(xyz)/6; 
% end 
% 
% 
% %% calculate volume integral 
% B = zeros(size(heart_coarse_surface.points,1)+size(torso_coarse_surface.points,1),size(NablaVm,3));
% ci = [0.1 0.2];
% co = [0.2 0];
% cin = 0.05;
% for surfaceID = 1:2
%     
%     if surfaceID == 1 % heart surface
%         Point_lookup = heart_coarse_surface.points;
%     else
%         Point_lookup = torso_coarse_surface.points;
%     end
%     
%     b = -(cin)/(2*pi*(ci(1,surfaceID)+co(1,surfaceID)));
%     R_pts = cellCenters_w_fib;
%     for pointID = 1:size(Point_lookup,1)
%         P_pts = repmat(Point_lookup(pointID,:), size(R_pts,1),1);
%         if surfaceID == 1
%             R_d = D(:,pointID);
%         else
%             R_d = D_torso(:,pointID);
%         end
% 
%         Diff_v = P_pts-R_pts;
%         
%         % calculate volume integral
%         for t_step = 1:size(NablaVm,3)
%             %pairwise multiplication.
%             NablaVm_t_i = squeeze(NablaVm_coarse(:,:,t_step));
%             if surfaceID == 1
%                 B(pointID,t_step) = b.*sum(sum(NablaVm_t_i.*Diff_v,2)./R_d.^3.*vol_elem,1);
%             else
%                 B(pointID+size(heart_coarse_surface.points,1),t_step) = b.*sum(sum(NablaVm_t_i.*Diff_v,2)./R_d.^3.*vol_elem,1);
%             end
%         end
%     end
% end
% 
% %% compute BSPM
% 
% heart_mesh = PrepareTriangleMesh(heart_coarse_surface.points,heart_coarse_surface.cells);
% torso_mesh = PrepareTriangleMesh(torso_coarse_surface.points,torso_coarse_surface.cells);
% meshes = {heart_mesh,torso_mesh};
% [Linv,~,~,~]=TransferMatrix_TMV_Linear(meshes, ci,co, cin, 2);
% 
% 
% bspm = Linv\B;
% 
% 
