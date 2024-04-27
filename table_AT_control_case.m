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
paths = {'/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1000_e8x8/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1000_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1000_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1500_e8x8/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1500_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r02_d1500_e16x16/2021-06-11_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r04_d1000_e8x8 /2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r04_d1000_e12x12/2021-06-11_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r04_d1000_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r04_d1500_e8x8/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r04_d1500_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r04_d1500_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r06_d1000_e8x8 /2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r06_d1000_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r06_d1000_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r06_d1500_e8x8/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r06_d1500_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r06_d1500_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r08_d1000_e8x8 /2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r08_d1000_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r08_d1000_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r08_d1500_e8x8/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r08_d1500_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r08_d1500_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r10_d1000_e8x8 /2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r10_d1000_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r10_d1000_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r10_d1500_e8x8/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r10_d1500_e12x12/2021-06-10_basic_1000.0/phie_recovery.igb';
    '/Volumes/Daten/Benutzer/pt732/Experiments/01_control_case/08_Simulations_parameters/r10_d1500_e16x16/2021-06-10_basic_1000.0/phie_recovery.igb'};

%% Matrix B


for i = 1:size(paths,1)

    phie = igb_data(paths{i}, 'all', 'all');
    
    B = cat(2,zeros(size(phie,1),5),phie);
    refSig = zeros(size(B,1),1);
    B = B - repmat(refSig, 1, size(B,2));
    
    B = phie; 
%     B = B(:,4:end);
    B = B(:,1:125);
    
    F = gradient(B');
    [min_F,ATF] = min(F);

    AT(i,1) = min(ATF);
    AT(i,2) = max(ATF);
end

varNames = {'AT_min';'AT_max'};

rowNames = {'r02_d1000_e08x08';'r02_d1000_e12x12';'r02_d1000_e16x16';'r02_d1500_e08x08';
    'r02_d1500_e12x12';'r02_d1500_e16x16';'r04_d1000_e08x08';'r04_d1000_e12x12';
    'r04_d1000_e16x16';'r04_d1500_e08x08';'r04_d1500_e12x12';'r04_d1500_e16x16';
    'r06_d1000_e08x08';'r06_d1000_e12x12';'r06_d1000_e16x16';'r06_d1500_e08x08';
    'r06_d1500_e12x12';'r06_d1500_e16x16';'r08_d1000_e08x08';'r08_d1000_e12x12';
    'r08_d1000_e16x16';'r08_d1500_e08x08';'r08_d1500_e12x12';'r08_d1500_e16x16';
    'r10_d1000_e08x08';'r10_d1000_e12x12';'r10_d1000_e16x16';'r10_d1500_e08x08';
    'r10_d1500_e12x12';'r10_d1500_e16x16'}; %1000 == 500 /// 1500 == 1000


T_AT = table(AT(:,1),AT(:,2),'VariableNames',varNames,'RowNames',rowNames);

% writetable(T_AT,'Table_AT.xlsx','Sheet','1','WriteVariableNames',true,'WriteRowNames',true);
