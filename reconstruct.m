clc
close all
clear

addpath(genpath('functions'));


gitpath = '/Volumes/bordeaux/IBT/matlab';
gitpath_m_p = strcat(gitpath,'/projects/');

gitpath = '/Volumes/bordeaux/IBT/matlab';
gitpath_m_c = strcat(gitpath,'/common/');

addpath(genpath(gitpath_m_c));

%%
load('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/intracardiac_inverse/elecTissDist2.1/matrices/A.mat')
load('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/intracardiac_inverse/elecTissDist2.1/matrices/L.mat')
load('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/intracardiac_inverse/elecTissDist2.1/signals/B.mat');

B = B(:,51:700);
B = B-repmat(mean(B,1),size(B,1),1);
A = A-repmat(mean(A,1),size(A,1),1);

% lambda = tikhonovLcurve_corner(B, A, L, true, [-6 3], 50, 1, gca);
lambda = 1e-3;

X = tikhonov(B, A, L, lambda, true);

%%
sigma = 20;
sig = gaussFiltfilt(X, sigma);
derivSig = (sig(:,3:end)-sig(:,1:end-2))/2;
[~,at] = max(derivSig,[],2);

%%
tiss = vtkRead('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/intracardiac_inverse/elecTissDist2.1/geometries/tissue.ply');
elec = vtkRead('/Volumes/Daten/Benutzer/pt732/Experiments/MATLAB/Codes/inverse/intracardiac_inverse/elecTissDist2.1/geometries/electrodes.ply');
comb = vtkAppendPolyData({tiss, elec});

qtrip(comb.points, comb.cells, X, 1041, [0 90 90], 10, [-80 80], 20);
pause(0.5)
qtrip(comb.points, comb.cells, at, 1042, [0 90 90], 20, [0 250], 20);
