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
paths = {'/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/8x8/tissue_ints_10_f1/2021-07-07_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/8x8/tissue_ints_20_f1/2021-07-07_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/8x8/tissue_ints_40_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/8x8/tissue_ints_60_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/12x12/tissue_ints_10_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/12x12/tissue_ints_20_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/12x12/tissue_ints_40_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/12x12/tissue_ints_60_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/16x16/tissue_ints_10_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/16x16/tissue_ints_20_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/16x16/tissue_ints_40_f1/2021-07-08_basic_2000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/05_fibrotic_reentry_case/simulations_d500/16x16/tissue_ints_60_f1/2021-07-08_basic_2000.0/phie_recovery.igb'};

%% Matrix B


for i = 1:size(paths,1)
    phie = igb_data(paths{i}, 'all', 'all');
    
    B = phie;  
    
    sampleFrequency = 333;
    limit = round(0.245/2*sampleFrequency);
    B = B(:,limit:end);
    
    F = gradient(B');
    [~,ATF] = max(F);

    AT(i,1) = min(ATF);
    AT(i,2) = max(ATF);
end

varNames = {'AT_min';'AT_max'};

rowNames = {'10_2_d500_e08x08';'20_2_d500_e08x08';'40_2_d500_e08x08';
    '60_2_d500_e08x08';'10_2_d500_e12x12';'20_2_d500_e12x12';
    '40_2_d500_e12x12';'60_2_d500_e12x12';'10_2_d500_e16x16';
    '20_2_d500_e16x16';'40_2_d500_e16x16';'60_2_d500_e16x16'}; 


T_AT = table(AT(:,1),AT(:,2),'VariableNames',varNames,'RowNames',rowNames);

% writetable(T_AT,'Table_AT.xlsx','Sheet','1','WriteVariableNames',true,'WriteRowNames',true);
