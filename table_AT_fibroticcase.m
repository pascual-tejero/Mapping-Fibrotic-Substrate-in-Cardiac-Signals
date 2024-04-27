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


%% Paths
paths = {'/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/8x8/tissue_ints_10_2/myiofibroblasts_2/2021-06-22_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/8x8/tissue_ints_20_2/2021-06-22_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/8x8/tissue_ints_40_2/2021-06-22_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/8x8/tissue_ints_60_2/2021-06-22_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/12x12/tissue_ints_10_2/2021-06-23_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/12x12/tissue_ints_20_2/2021-06-23_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/12x12/tissue_ints_40_2/2021-06-23_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/12x12/tissue_ints_60_2/2021-06-23_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/16x16/tissue_ints_10_2/2021-06-28_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/16x16/tissue_ints_20_2/2021-06-28_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/16x16/tissue_ints_40_2/2021-06-28_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/04_fibrotic_tissue_v2/16x16/tissue_ints_60_2/2021-06-28_basic_1000.0/phie_recovery.igb'};

%% Matrix B


for i = 1:size(paths,1)
    phie = igb_data(paths{i}, 'all', 'all');
    
    B = cat(2,zeros(size(phie,1),5),phie);
    refSig = zeros(size(B,1),1);
    B = B - repmat(refSig, 1, size(B,2));
    
    B = B(:,4:end);
    B = B(:,1:140);
    
    F = gradient(B');
    [min_F,ATF] = min(F);

    AT(i,1) = min(ATF);
    AT(i,2) = max(ATF);
end

varNames = {'AT_min';'AT_max'};

rowNames = {'10_2_d500_e08x08';'20_2_d500_e08x08';'40_2_d500_e08x08';
    '60_2_d500_e08x08';'10_2_d500_e12x12';'20_2_d500_e12x12';
    '40_2_d500_e12x12';'60_2_d500_e12x12';'10_2_d500_e16x16';
    '20_2_d500_e16x16';'40_2_d500_e16x16';'60_2_d500_e16x16'};


T_AT = table(AT(:,1),AT(:,2),'VariableNames',varNames,'RowNames',rowNames);
% 
% writetable(T_AT,'Table_AT.xlsx','Sheet','1','WriteVariableNames',true,'WriteRowNames',true);

