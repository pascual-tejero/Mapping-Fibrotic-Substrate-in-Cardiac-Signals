clc
clear
close all

%% Add folders to search path that are necessary for running the simulation 

addpath('/Volumes/bordeaux/IBT10.10/matlab/common/vtkToolbox/MATLAB')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/inverse_problem/BEM')
addpath('/Volumes/bordeaux/IBT10.10/matlab/thirdparty/gptoolbox')
addpath('/Volumes/beaujolais/Benutzer/js191/Pascual/tikhonov')
addpath('/Volumes/beaujolais/Benutzer/js191/Pascual/operators')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/EGMDataProcessing/FilteringAndProcessing')
addpath('/Volumes/bordeaux/IBT10.10/matlab/projects/EGMDataProcessing/activityAndDPs')
%addpath('activationtimes-master/computeGradientAndHessian')

%% Tissue mesh

xo_points=0:1000:50000; % x axis
yo_points=0:1000:50000; % y axis
zo_points=800:200:1800; % z axis

% Creation of a 3-D meshgrid
[x_points, y_points, z_points] = meshgrid(xo_points,yo_points,zo_points);

% Creation of a bounding area or volume that envelops a set of 3-D points
shp_points = alphaShape(x_points(:),y_points(:),z_points(:)); 
% plot(shp_points)

% Creation of a matrix representing the facets that make up the boundary of the alpha shape
[bf_points,P_points] = boundaryFacets(shp_points); 

% Representation of any 2-D or 3-D triangulation data that is in matrix format
TR_points = triangulation(bf_points,P_points);

% Exportation of triangulations and triangulated surfaces to vtk file format.
vtk_coarse = triangulationToVtk(TR_points);

% Saving the vtk file
vtkWrite(vtk_coarse,'sources.vtk',false,'ascii');

%% Electrodes mesh

% Creation first of the position where the electrode will be
distance_electrodes = 3800; 

% 8x8 electrode distribution
% xo_elec=14500:3000:35500;  
% yo_elec=14500:3000:35500;

% 12x12 electrode distribution
xo_elec=8500:3000:41500; 
yo_elec=8500:3000:41500;

% 16x16 electrode distribution
% xo_elec=2500:3000:47500;
% yo_elec=2500:3000:47500;

% cont = 0;
% for i=1:length(xo_elec)
%     electrodes_distribution(i,:) = [(1+cont):(length(xo_elec)+cont)];
%     cont = length(xo_elec)*i;
% end

% Number of electrodes used in this simulation
number_electrodes = length(xo_elec)*length(yo_elec);

% Creation of a 3-D meshgrid
[x_elec, y_elec] = meshgrid(xo_elec,yo_elec); %Return a 3D grid, according to the inputs

% Changing the z-axis distance of the electrodes so that they are exactly 0.5 mm from the tissue mesh

points_electrodes = [x_elec(:),y_elec(:),distance_electrodes.*ones(length(xo_elec)*length(yo_elec),1)];
 
% Creation of 4-point spheres
[x_sphere,y_sphere,z_sphere] = sphere(4);

r = 1; %Radio of spheres

X_sphere = x_sphere * r;
Y_sphere = y_sphere * r;
Z_sphere = z_sphere * r;

X_sphere = X_sphere(:);
Y_sphere = Y_sphere(:);
Z_sphere = Z_sphere(:);

electrode_sphere = []; 

for i=1:size(points_electrodes,1)
    X_sphere2 = points_electrodes(i,1) + X_sphere;
    Y_sphere2 = points_electrodes(i,2) + Y_sphere;
    Z_sphere2 = points_electrodes(i,3) + Z_sphere;

    temp = [X_sphere2 Y_sphere2 Z_sphere2];
    electrode_sphere = [electrode_sphere; temp];
end

% Creation of a bounding area or volume that envelops a set of 3-D points
sphere_elec = alphaShape(electrode_sphere(:,1),electrode_sphere(:,2),electrode_sphere(:,3));
% plot(shpere_elec)

% Creation of a matrix representing the facets that make up the boundary of the alpha shape
[bf_elec,P_elec] = boundaryFacets(sphere_elec); 

% Representation of any 2-D or 3-D triangulation data that is in matrix format
TR_elec = triangulation(bf_elec,P_elec);

% Exportation of triangulations and triangulated surfaces to vtk file format.
vtk_electrodes = triangulationToVtk(TR_elec); 

% Saving the vtk file
vtkWrite(vtk_electrodes,'electrodes.vtk',false,'ascii');


%% Bounding box

% Limits of bounding box
xo_bb=-20000:10000:80000;
yo_bb=-20000:10000:80000;
zo_bb=600:1000:2600;

% Creation of a 3-D meshgrid
[x_bb, y_bb, z_bb] = meshgrid(xo_bb,yo_bb,zo_bb); 

% Creation of a bounding area or volume that envelops a set of 3-D points
shp_bb = alphaShape(x_bb(:),y_bb(:),z_bb(:)); %Creates a bounding area or volume that envelops a set of 3-D points

% Creation of a matrix representing the facets that make up the boundary of the alpha shape
[bf_bb,P_bb] = boundaryFacets(shp_bb);

% Representation of any 2-D or 3-D triangulation data that is in matrix format
TR_bb = triangulation(bf_bb,P_bb);

% Exportation of triangulations and triangulated surfaces to vtk file format.
vtk_coarse_bbox = triangulationToVtk(TR_bb);


%% Transfer matrix
% 
% if exist('A') % If exists tranfer matrix A...
%     load('A.mat');
% else % Otherwise, it is created with the function below
    A = computeTransferMat_A2(vtk_coarse, vtk_electrodes,vtk_coarse_bbox);  
% end
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/A.mat', 'A', '-v7.3');

%% Matrix L (regularization term)

% Calculation of matrix L by using cotangent-Laplacian method
L = Laplacian_transmural(vtk_coarse);

% Calculation of matrix L by using weighting of the regularization of cotangent-Laplacian method
% A: Lead field matrix
% 
% L: Laplacian matrix
% 

p = 0.5; % try different values

w = vecnorm(A, 2, 1)';

w = (w/max(w)).^p;

W = spdiags(w, 0, size(A,2), size(A,2));

L_weighted = W*L;

% Calculation of matrix L by using laplacian based on "locally quadratic surface fits"
% geom = vtk_coarse;
% geom.node = vtk_coarse.points;
% geom.face = vtk_coarse.cells;
% 
% w = wghFcn(1);
% 
% [cDf cHf] = meshVolDiffHessMatrix(geom,wghFcn)

% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/L.mat', 'L', '-v7.3');

%% Matrix B (electrograms)

sampleFrequency = 334;

%Load electrogram data
phie = igb_data('/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations/12x12/tissue_ints_60_f1/dis_source-electrodes_10000mm /2021-09-14_basic_2000.0/phie_recovery.igb', 'all', 'all');
% phie = [phie phie phie phie];

B = phie;
% B = movmean(B,[6 6],2); %Filtering the reconstructed signal with

%Cut the signal until the second stimulus (~240ms)
% limit = round(0.245/2*sampleFrequency);
limit = 200;
cut_simulation(1) = 2000;
cut_simulation(2) = round(limit/(2*sampleFrequency)*2000);
B = B(:,limit:end);

% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/B.mat', 'B', '-v7.3');

%% Calculation of Activation Time (AT) and Reporalization Time (RT) of the values of the electrograms (ground truth)

% Function to calculate AT and RT
[AT_groundtruth,RT_groundtruth] = calculation_AT_RT(B,cut_simulation,20);

% Ordering of the AT and RT matrix
cont = 0;
for i=1:(length(AT_groundtruth)/length(xo_elec))
    values = i:length(xo_elec):(length(AT_groundtruth)-length(xo_elec)+i);
    
    AT_groundtruth_new(1+cont:length(xo_elec)+cont) = AT_groundtruth(values);
    RT_groundtruth_new(1+cont:length(xo_elec)+cont) = RT_groundtruth(values);
    
    cont = length(xo_elec)*i;
end

%% Inverse problem using Tikhonov regularization

% Function to use Tikhonov regularization
[U,Q,S,M] = tikhonov_gsvd_decompose(A,L_weighted);
% [U,Q,S,M] = tikhonov_gsvd_decompose(A,L_weighted);

% Selection of the correct regularization term
lambda = tikhonov_gsvd_Lcurve(B,U,Q,S,M,[-5 2]);
lambda = lambda/10;
[X,~,~] = tikhonov_gsvd(B,U,Q,S,M,lambda);
% [X,resNorm,solNorm] = tikhonov(B, A, L, 1e-5);

%Xm = movmean(X,[6 6],2); %Filtering the reconstructed signal with

disp('Smoothing...')
sigma = 1.1;
Xm = gaussFiltfilt(X, sigma); % Gauss filter
% derivSig = gradient(sig); %derivSig = (sig(:,3:end)-sig(:,1:end-2))/2; %
% [~,at] = max((derivSig),[],2);

%% Calculation AT and RT - Reconstructed matrix Xm
[AT_reconstruction,RT_reconstruction] = calculation_AT_RT(Xm,cut_simulation,20);
figure; plot(AT_reconstruction, '+'); hold on; plot(RT_reconstruction, '+');

y_lins = linspace(cut_simulation(2),cut_simulation(1),size(Xm,2));
% figure; plot(y_lins,B(1,:)); hold on; plot(y_lins,Xm(20,:));

%% Calculation of Root Mean Square Error (RMSE) of AT and RT (between groundtruth and reconstruction)

% Find the closest points to the electrodes
distance = pdist2(vtk_coarse.points, points_electrodes);
[distance_min, distance_pos_min] = min(distance);

AT_reconstruction_pos = AT_reconstruction(distance_pos_min); AT_reconstruction_pos2 = AT_reconstruction_pos;
RT_reconstruction_pos = RT_reconstruction(distance_pos_min); RT_reconstruction_pos2 = RT_reconstruction_pos;


mean_AT_reconstruction_pos = mean(AT_groundtruth);
std_AT_reconstruction_pos = std(AT_groundtruth);

mean_RT_reconstruction_pos = mean(RT_groundtruth);
std_RT_reconstruction_pos = std(RT_groundtruth);


AT_reconstruction_pos2((AT_reconstruction_pos>=mean_AT_reconstruction_pos+std_AT_reconstruction_pos) | (AT_reconstruction_pos<=mean_AT_reconstruction_pos-std_AT_reconstruction_pos)) = NaN;
RT_reconstruction_pos2((RT_reconstruction_pos>=mean_RT_reconstruction_pos+std_RT_reconstruction_pos) | (RT_reconstruction_pos<=mean_RT_reconstruction_pos-std_RT_reconstruction_pos)) = NaN;

val = 1:4:length(xo_elec);

cont1 = 1:length(xo_elec);
cont2 = 1;

% RMSE
RMSE_AT = sqrt((AT_reconstruction_pos2 - AT_groundtruth).^2);
RMSE_RT = sqrt((RT_reconstruction_pos2 - RT_groundtruth).^2);

% Reshape of RMSE
RMSE_reshape_AT = reshape(RMSE_AT,length(xo_elec),length(yo_elec)); RMSE_reshape_AT_n = RMSE_reshape_AT';
RMSE_reshape_RT = reshape(RMSE_RT,length(xo_elec),length(yo_elec)); RMSE_reshape_RT_n = RMSE_reshape_RT';

for i = 1:length(xo_elec)/4
    for j = 1:length(yo_elec)/4
        
        values_selected_AT(val(i):(val(i)+3),val(j):(val(j)+3)) = RMSE_reshape_AT_n(val(i):(val(i)+3),val(j):(val(j)+3));
        values_selected_RT(val(i):(val(i)+3),val(j):(val(j)+3)) = RMSE_reshape_RT_n(val(i):(val(i)+3),val(j):(val(j)+3));
        
        % Mean RMSE
        RMSE_mean_AT(cont1(cont2)) = mean(values_selected_AT(val(i):(val(i)+3),val(j):(val(j)+3)),'all','omitnan');
        RMSE_mean_RT(cont1(cont2)) = mean(values_selected_RT(val(i):(val(i)+3),val(j):(val(j)+3)),'all','omitnan');
       
        % Standard deviation RMSE
        RMSE_std_AT(cont1(cont2)) = nanstd(values_selected_AT(val(i):(val(i)+3),val(j):(val(j)+3)),1,'all');
        RMSE_std_RT(cont1(cont2)) = nanstd(values_selected_RT(val(i):(val(i)+3),val(j):(val(j)+3)),1,'all');
         
       cont2 = cont2 + 1;

    end
end
%%
filename = 'case_reentry_fibrotic_d10_12x12_60col.xlsx';

writematrix(RMSE_mean_AT,filename,'Sheet','RMSE_mean_AT');

writematrix(RMSE_mean_RT,filename,'Sheet','RMSE_mean_RT');

writematrix(RMSE_std_AT,filename,'Sheet','RMSE_std_AT');

writematrix(RMSE_std_RT,filename,'Sheet','RMSE_std_RT');

%% Representation of reconstruction and, AT and RT map
comb = vtkAppendPolyData({vtk_coarse, vtk_electrodes});
%%
% surf=alphaShape(points_electrodes(:,1),points_electrodes(:,2)); %2000.*ones(length(y_points(:)),1))
% [trisurf,psurf]=alphaTriangulation(surf);
% psurf(:,3)=points_electrodes(:,3);

% qtrip(comb.points, comb.cells, Xm, 1041, [0 90 90], 2, [min(min(Xm)) max(max(Xm))], 10); %Reconstruccion potencial extracelular
% pause(0.5)
% qtrip(comb.points, comb.cells, AT_reconstruction', 1042, [0 90 90], 10, [min(min(AT_reconstruction)) max(max(AT_reconstruction))]);
% pause(0.5)
% qtrip(comb.points, comb.cells, RT_reconstruction', 1042, [0 90 90], 10, [min(min(RT_reconstruction)) max(max(RT_reconstruction))]);

% pause(0.5)
% qtrip(psurf, trisurf, AT_groundtruth_new', 1042, [0 90 90], 5, [min(min(AT_groundtruth_new)) max(max(AT_groundtruth_new))]);
% pause(0.5)
% qtrip(psurf, trisurf, RT_groundtruth_new', 1042, [0 90 90], 2, [min(min(RT_groundtruth)) max(max(RT_groundtruth))]);
% 
% figure;
% imagesc(reshape(AT_groundtruth_new,length(xo_elec),length(yo_elec)))
% figure
% imagesc(reshape(RT_groundtruth_new,length(xo_elec),length(yo_elec)))
% 
%%
[singularity_points,phase] = singularity_points_function(Xm, vtk_coarse);

qtrip(comb.points, comb.cells, phase, 1041, [0 90 90], 0.5, [min(min(phase)) max(max(phase))], 10);
qtrip(comb.points, comb.cells, singularity_points, 1041, [0 90 90], 0.5, [0 1], 10);

% fibroti%%
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/AT_reconstruction.mat', 'AT_reconstruction', '-v7.3');
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/vtk_coarse.mat', 'vtk_coarse', '-v7.3');

% Singularidad 16x16 --> 10% fibrotico y 60% fibrotico, reentrada y tambien tejido sano 16x16 (puede que 8x8 para ver que no es muy optimo)
% Distancia 2 mm --> 16x16 electrodos, 10% y 60% fibrotico

%% Reconstruction with Tikhonov-2nd-Order:
% lambda=0.001;
% x_reconstruction = (A'*A+lambda^2*(L_weighted'*L_weighted))\ A*Xm;
% 
% m_x_reconstruction = mean(x_reconstruction,1);
% for i=1:size(x_reconstruction,1)
% x_reconstruction(i,:) = x_reconstruction(i,:)-m_x_reconstruction; % remove mean over points in heart, time-step-wise
% end
% 
% save x_reconstruction.dat x_reconstruction -ascii

