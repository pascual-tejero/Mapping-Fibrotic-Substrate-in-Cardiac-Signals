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
addpath('activationtimes-master/computeGradientAndHessian')


%% Tissue mesh

xo_points=-25000:1000:25000; % x axis
yo_points=-25000:1000:25000; % y axis
zo_points=-500:1000:500; % z axis

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
distance_electrodes = 3500;
% 
xo_elec=-10500:3000:10500;
yo_elec=-10500:3000:10500;
% % 
% xo_elec=-16500:3000:16500;
% yo_elec=-16500:3000:16500;
% % 
% xo_elec=-22500:3000:22500;
% yo_elec=-22500:3000:22500;

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
 %1000 or 1500
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
shpere_elec = alphaShape(electrode_sphere(:,1),electrode_sphere(:,2),electrode_sphere(:,3));
% plot(shpere_elec)

% Creation of a matrix representing the facets that make up the boundary of the alpha shape
[bf_elec,P_elec] = boundaryFacets(shpere_elec); 

% Representation of any 2-D or 3-D triangulation data that is in matrix format
TR_elec = triangulation(bf_elec,P_elec);

% Exportation of triangulations and triangulated surfaces to vtk file format.
vtk_electrodes = triangulationToVtk(TR_elec); 

% Saving the vtk file
vtkWrite(vtk_electrodes,'electrodes.vtk',false,'ascii');


%% Bounding box

% Limits of bounding box
xo_bb=-30000:5000:30000;
yo_bb=-30000:5000:30000;
zo_bb=-1000:1000:2000;

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

A = computeTransferMat_A2(vtk_coarse, vtk_electrodes,vtk_coarse_bbox);  
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/A.mat', 'A', '-v7.3');

%% Matrix L (regularization term)

% Calculation of matrix L by using cotangent-Laplacian method
L = Laplacian_transmural(vtk_coarse);

% Calculation of matrix L by using weighting of the regularization of cotangent-Laplacian method
% A: Lead field matrix
% 
% L: Laplacian matrix
% 
p = 10; % try different values

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

%% Matrix B
sampleFrequency = 334;

%Load electrogram data
phie = igb_data('/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations/8x8/tissue_ints_10_f1/dis_source-electrodes_2mm/2021-09-14_basic_2000.0/phie_recovery.igb', 'all', 'all');

B_reconstruction = phie;
B_reconstruction = B_reconstruction(:,1:150);

B = [B_reconstruction B_reconstruction B_reconstruction];


%Cut the signal 
limit = 150;
cut_simulation(1) = 3000;
cut_simulation(2) = round(limit/(3*sampleFrequency)*3000);
% 
B = B(:,limit:end);

%% Calculation of Activation Time (AT) and Reporalization Time (RT) of the values of the electrograms (ground truth)

% Function to calculate AT and RT
[AT_groundtruth,RT_groundtruth] = calculation_AT_RT(B,cut_simulation,6);

% Ordering of the AT and RT matrix
cont = 0;
for i=1:(length(AT_groundtruth)/length(xo_elec))
    values = i:length(xo_elec):(length(AT_groundtruth)-length(xo_elec)+i);
    
    AT_groundtruth_new(1+cont:length(xo_elec)+cont) = AT_groundtruth(values);
    RT_groundtruth_new(1+cont:length(xo_elec)+cont) = RT_groundtruth(values);
    
    cont = length(xo_elec)*i;
end
% figure; plot(AT_groundtruth_new, '+'); hold on; plot(RT_groundtruth_new, '+');

figure; plot(AT_groundtruth, '+'); hold on; plot(RT_groundtruth, '+');

%% Inverse problem using Tikhonov regularization

% Function to use Tikhonov regularization
[U,Q,S,M] = tikhonov_gsvd_decompose(A,L);

% Selection of the correct regularization term
lambda = tikhonov_gsvd_Lcurve(B_reconstruction,U,Q,S,M,[-5 2]);
[X,~,~] = tikhonov_gsvd(B_reconstruction,U,Q,S,M,lambda);

Xm = movmean(X,[3 3],2);%Filtering the reconstructed signal with

%% Calculation AT and RT - Reconstructed matrix Xm
Xm_AT_RT = [Xm Xm Xm];

Xm_AT_RT = Xm_AT_RT(:,limit:end);

[AT_reconstruction,RT_reconstruction] = calculation_AT_RT(Xm_AT_RT,cut_simulation,20);

figure; plot(AT_reconstruction, '+'); hold on; plot(RT_reconstruction, '+');

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
% filename = 'case_healthytissue_reentry_d0.5_12x12.xlsx';
% 
% writematrix(RMSE_mean_AT,filename,'Sheet','RMSE_mean_AT');
% 
% writematrix(RMSE_mean_RT,filename,'Sheet','RMSE_mean_RT');
% 
% writematrix(RMSE_std_AT,filename,'Sheet','RMSE_std_AT');
% 
% writematrix(RMSE_std_RT,filename,'Sheet','RMSE_std_RT');

%% Representation of reconstruction and, AT and RT map
comb = vtkAppendPolyData({vtk_coarse, vtk_electrodes});
% 
% surf=alphaShape(points_electrodes(:,1),points_electrodes(:,2)); %2000.*ones(length(y_points(:)),1))
% [trisurf,psurf]=alphaTriangulation(surf);
% psurf(:,3)=points_electrodes(:,3);
% 
% 
qtrip(comb.points, comb.cells, Xm, 1041, [0 90 90], 0.5, [round(min(min(Xm))) round(max(max(Xm)))], 15);
pause(0.5)
qtrip(comb.points, comb.cells, AT_reconstruction', 1042, [0 90 90], 20, [round(min(min(AT_reconstruction))) round(max(max(AT_reconstruction)))]);
pause(0.5)
qtrip(comb.points, comb.cells, RT_reconstruction', 1042, [0 90 90], 20, [round(min(min(RT_reconstruction))) round(max(max(RT_reconstruction)))]);

% pause(0.5)
% qtrip(psurf, trisurf, AT_groundtruth_new', 1042, [0 90 90], 5, [min(min(AT_groundtruth_new)) max(max(AT_groundtruth_new))]);
% pause(0.5)
% qtrip(psurf, trisurf, RT_groundtruth_new', 1042, [0 90 90], 2, [min(min(RT_groundtruth)) max(max(RT_groundtruth))]);
% figure
% imagesc(reshape(AT_groundtruth_new,length(xo_elec),length(yo_elec)))
% figure
% imagesc(reshape(RT_groundtruth_new,length(xo_elec),length(yo_elec)))