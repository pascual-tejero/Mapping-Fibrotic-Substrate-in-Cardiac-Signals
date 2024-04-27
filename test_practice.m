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

%% Initial points (mesh)
xo=0:1000:50000;
yo=0:1000:50000;
zo=800:400:2000;

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

vtkWrite(vtk_coarse,'sources.vtk',false,'ascii');

%% Electrodes mesh

xoe=14500:3000:35500;
yoe=14500:3000:35500;

% xoe=8500:3000:41500;
% yoe=8500:3000:41500;

% xoe=2500:3000:47500;
% yoe=2500:3000:47500;

cont = 0;
for i=1:length(xoe)
    electrodes_distribution(i,:) = [(1+cont):(length(xoe)+cont)];
    cont = length(xoe)*i;
end


[x, y] = meshgrid(xoe,yoe); %Return a 3D grid, according to the inputs

points_electrodes=[x(:),y(:),2500.*ones(length(xoe)*length(yoe),1)]; %distance = 1000 or 1500 
 
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

for i=1:size(points_electrodes,1)
    X_sphere = points_electrodes(i,1) + Xr;
    Y_sphere = points_electrodes(i,2) + Yr;
    Z_sphere = points_electrodes(i,3) + Zr;

    temp = [X_sphere Y_sphere Z_sphere];
    electrode_sphere = [electrode_sphere; temp];
end

shp_elec = alphaShape(electrode_sphere(:,1),electrode_sphere(:,2),electrode_sphere(:,3));
% plot(shp_elec);

% boundaryFacets
%Boundary facets bf, returned as a matrix. bf is of size m-by-n, where m is 
%the number of boundary facets and n is the number of vertices per facet.
%Vertex coordinates P, returned as a matrix. P is of size N-by-dim, where N 
%is the number of points on the boundary of the alpha shape and dim is either
%2 or 3 (for either a 2-D or 3-D alpha shape).

[bf,P] = boundaryFacets(shp_elec); %Returns a matrix representing the facets that 
%make up the boundary of the alpha shape

% triangulation
%Use triangulation to create an in-memory representation of any 2-D or 3-D 
%triangulation data that is in matrix format
TR_elec = triangulation(bf,P);

vtk_electrodes = triangulationToVtk(TR_elec); %vtkfile

vtkWrite(vtk_electrodes,'electrodes.vtk',false,'ascii');


%% Bounding box
xob=-30000:5000:30000;
yob=-30000:5000:30000;
zob=-1000:1000:2000;

[xb, yb, zb] = meshgrid(xob,yob,zob); %Return a 3D grid, according to the inputs

shpb = alphaShape(xb(:),yb(:),zb(:)); %Creates a bounding area or volume that envelops a set of 3-D points

[bfb,Pb] = boundaryFacets(shpb); %Returns a matrix representing the facets that 

TRb = triangulation(bfb,Pb);
vtk_coarse_bbox = triangulationToVtk(TRb);


%% Transfer matrix
exist_A = exist('A');

if exist_A == 0
    A = computeTransferMat_A2(vtk_coarse, vtk_electrodes,vtk_coarse_bbox);  %This needs to be calculated before hand TODO create a function to generate the A matrix if not exists
else
    load('A.mat');
end
% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/A.mat', 'A', '-v7.3');

%% Matrix L (regularization term)
L = Laplacian_transmural(vtk_coarse);

% A: Lead field matrix

% L: Laplacian matrix

% p = 2; % try different values
% 
% w = vecnorm(A, 2, 1)';
% 
% w = (w/max(w)).^p;
% 
% W = spdiags(w, 0, size(A,2), size(A,2));
% 
% L_weighted = W*L;

% save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/L.mat', 'L', '-v7.3');
%% Matrix B
sampleFrequency = 333;
%Load electrogram data
phie = igb_data('/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/8x8/tissue_ints_10_f1/fibrotic_tissue_reentry/2021-07-16_basic_2000.0/phie_recovery.igb', 'all', 'all');

B = phie;

%Cut the signal until the second stimulus (~240ms)
limit = round(0.245/2*sampleFrequency);

B = B(:,limit:end);

save('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/MATLAB/B.mat', 'B', '-v7.3');

%% Selection of a time window - Groundtruth
matrix = B;

sampleFrequency = 333;
mode = 1;
k = 1;
lenghtOfWin = 1;

%Obtention of the time window
accp = zeros(627,64);

for i = 1:size(matrix,1)
%     [stepFuncFinal(:,i),segments{:,i}] = getActiveSegmentsFromNLEOuni(signalOut(:,i),sampleFrequency,lenghtOfWin,k);
    [accp,timelngth] = Hactseg(matrix(i,:),sampleFrequency);
    accpt(i,1:length(accp)) = accp;
end

% Obtaining the AT of the last window
%Location of the last time window in order to calculate the AT
y_lins = linspace(40,2000,size(matrix,2));
window_s = zeros(1,20);

for i=1:size(accpt,1)
    pos_last1_accpt = find(accpt(i,:),1,'last');
    vector_values = accpt(i,1:pos_last1_accpt);
    pos_first1_accpt = find(vector_values == 0);

    pos_window(i,1) = pos_first1_accpt(end);
    pos_window(i,2) = pos_last1_accpt;

    window_B = matrix(i,pos_window(i,1):pos_window(i,2));

    pos_window_s(i,1) = pos_window(i,1)/size(matrix,2)*2000;
    pos_window_s(i,2) = pos_window(i,2)/size(matrix,2)*2000;

    value_window_s = y_lins(pos_window(i,1):pos_window(i,2));
    window_s(i,1:length(value_window_s)) = value_window_s;

%         figure; plot(y_lins,matrix(i,:)); hold on; plot(y_lins,accpt(i,:)); hold on; 
%         plot(value_window_s, window_B); 
%         legend("Electrogram","Detection function","Last time window");
%         title(['Fibrotic tissue - Reentry case - Electrode: ', num2str(i)]); 
%         xlabel("Time (ms)"); ylabel("Amplitude voltage (mV)")

    derivative = gradient(window_B);
    [~,at(i)] = min((derivative),[],2);

    matrix_AT(i) = window_s(i,at(i));
end

 AT_reconstruction = matrix_AT;



%% Calculation of the reporalization times


% for i=1:size(accpt,1)
%     pos_last1_accpt = find(accpt(i,:),1,'last');
%     vector_values = accpt(i,1:pos_last1_accpt);
%     
% end

%% Calculation of RMSE (AT_Groundtruth and AT_Reconstruction)
distance = pdist2(vtk_coarse.points, points_electrodes);
[distance_min, distance_pos_min] = min(distance');

AT_reconstruction_pos = AT_reconstruction(distance_pos_min);


%RMSE 
RMSE = sqrt((AT_reconstruction_pos - AT_groundtruth).^2);        
RMSE_reshape = reshape(RMSE,length(xoe),length(yoe));
RMSE_reshape = RMSE_reshape';

val = 1:2:length(xoe);
for i = 1:length(xoe)/2
    for j = 1:length(yoe)/2
        values_selected = RMSE_reshape(val(i):(val(i)+1),val(j):(val(j)+1));
        RMSE_mean(val(i):(val(i)+1),val(j):(val(j)+1)) = mean(RMSE_reshape(val(i):(val(i)+1),val(j):(val(j)+1)),'all');
    end
end




% coldist = 2*ones(1,length(xoe)*length(yoe)/4);
% mat2cell(electrodes_distribution, 2, coldist);




