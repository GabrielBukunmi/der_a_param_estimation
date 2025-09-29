%% Generate paper plots and additional figures


% --- Case1/plot_EKF_UKF_Case1.m ------------------------------------------
thisDir = fileparts(mfilename('fullpath'));  
cd(thisDir);
disp('Running Case1/plot_EKF_UKF_Case1.m ...');
run(fullfile(thisDir, 'Case1', 'plot_EKF_UKF_Case1.m'));

% --- Case2/plot_EKF_UKF_Case2.m --------------
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
disp('Running Case2/plot_EKF_UKF_Case2.m ...');
run(fullfile(thisDir, 'Case2', 'plot_EKF_UKF_Case2.m'));

% --- der_a.m --------------------------------------------------------------
thisDir = fileparts(mfilename('fullpath'));  
cd(thisDir);                                 
disp('Running der_a.m ...');
run(fullfile(thisDir, 'der_a.m'));

%--- observability plots---------------------------------------------------
thisDir = fileparts(mfilename('fullpath'));  
cd(thisDir);                                 
disp('Running Observabilityplots ...');
run(fullfile(thisDir, 'Observability', 'plotsingularmatrix.m'));
run(fullfile(thisDir, 'Observability', 'param_weight_contributions.m'));

%---- Plot smooth saturation and smooth deadband functions-----------------
thisDir = fileparts(mfilename('fullpath'));  
cd(thisDir);                                 
disp('plotting saturation functions ...');
run(fullfile(thisDir, 'PLOTSSF.m'));

disp(' Done. All plots generated!!.');
