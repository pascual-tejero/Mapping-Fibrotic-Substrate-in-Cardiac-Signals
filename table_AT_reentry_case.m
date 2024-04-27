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
paths = {'/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02/d1000_e8x8/2021-08-03_basic_3000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02/d1500_e8x8/2021-08-03_basic_3000.0/phie_recovery.igb';   
    '/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02/d1000_e12x12/2021-08-03_basic_3000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02/d1500_e12x12/2021-08-03_basic_3000.0/phie_recovery.igb';  
    '/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02/d1000_e16x16/2021-08-03_basic_3000.0/phie_recovery.igb';    
    '/Volumes/Daten/Benutzer/pt732/Experiments/03_control_reentry/Simulations_r02/d1500_e16x16/2021-08-03_basic_3000.0/phie_recovery.igb'};

%% Matrix B


for i = 1:size(paths,1)
    phie = igb_data(paths{i}, 'all', 'all');
    
    B = phie;
    
    refSig = zeros(size(B,1),1);
    B = B - repmat(refSig, 1, size(B,2));
    
    B = B(:,end-300:end);
    
    F = gradient(B');
    [min_F,ATF] = min(F);

    AT(i,1) = min(ATF);
    AT(i,2) = max(ATF);
end

varNames = {'AT_min';'AT_max'};

rowNames = {'r02_d1000_e08x08';'r02_d1500_e08x08';'r02_d1000_e12x12';'r02_d1500_e12x12';
    'r02_d1000_e16x16';'r02_d1500_e16x16'}; %1000 == 500 /// 1500 == 1000


T_AT = table(AT(:,1),AT(:,2),'VariableNames',varNames,'RowNames',rowNames);

% writetable(T_AT,'Table_AT.xlsx','Sheet','1','WriteVariableNames',true,'WriteRowNames',true);

