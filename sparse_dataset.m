%% Random seed
seed = 1234;
rng(seed);

%% Data
CAMELS_US = readtable('data/CAMELS_US.csv');
[groupidx, group] = findgroups(CAMELS_US.catchment_id);

% Store catchment data in a cell, cell row number = catchment index
catchment_datas = cell(length(group),1);
for i = 1:length(group)
    catchment_datas{i,1} = CAMELS_US(groupidx==i, :);
end

%% model list
model_classes  = {
    'm_01_collie1_1p_1s',...
    'm_02_wetland_4p_1s',...
    'm_03_collie2_4p_1s',...
    'm_04_newzealand1_6p_1s',...
    'm_05_ihacres_7p_1s',...
    'm_06_alpine1_4p_2s',...
    'm_07_gr4j_4p_2s',...
    'm_08_us1_5p_2s',...
    'm_09_susannah1_6p_2s',...
    'm_10_susannah2_6p_2s',...
    'm_11_collie3_6p_2s',...
    'm_12_alpine2_6p_2s',...
    'm_13_hillslope_7p_2s',...
    'm_14_topmodel_7p_2s',...
    'm_15_plateau_8p_2s',...
    'm_16_newzealand2_8p_2s',...
    'm_17_penman_4p_3s',...
    'm_18_simhyd_7p_3s',...
    'm_19_australia_8p_3s',...
    'm_20_gsfb_8p_3s',...
    'm_21_flexb_9p_3s',...
    'm_22_vic_10p_3s',...
    'm_24_mopex1_5p_4s',...
    'm_25_tcm_6p_4s',...
    'm_26_flexi_10p_4s',...
    'm_27_tank_12p_4s',...
    'm_28_xinanjiang_12p_4s',...
    'm_29_hymod_5p_5s',...
    'm_30_mopex2_7p_5s',...
    'm_31_mopex3_8p_5s',...
    'm_32_mopex4_10p_5s',...
    'm_34_flexis_12p_5s',...
    'm_35_mopex5_12p_5s',...
    'm_36_modhydrolog_15p_5s',...
    'm_37_hbv_15p_5s'};

%% Create input grid for evaluation
n_model_class = length(model_classes);                       % number of models
n_para_sets = 1000;                                          % number of parameters sets
n_catchments = length(group);                                % number of catchments; max = 533
n_catchment_eval = round(n_catchments*0.05);                 % number of catchments evaluated for each model

n_eval = n_model_class * n_para_sets * n_catchment_eval;     % number of model evalution
eval_grid = cell(n_eval,4);                                  % creating 'eval_grid' to store input parameters

% iterate over model class; select random parameters for each model class
% to create model instance; select random catchments for each instance
for i = 1:length(model_classes)
    
    model_class = model_classes{i};                          % model_class of the iteration

    % Sampling model parameter
    % extract the parameter range for this model class
    model_range = feval([model_class,'_parameter_ranges']);  % Call the function that stores parameter ranges
    
    % Find number parameters
    numPar = size(model_range,1);                            % Identify the number of model parameters
    
    % Sample a parameter set from the range
    input_thetas = model_range(:,1)+rand(numPar,n_para_sets).*(model_range(:,2)-model_range(:,1)); % random parameters, numPar * n_para_sets
    input_thetas = transpose(input_thetas); % convert to shape of n_para_sets*numPar
    input_thetas = kron(input_thetas, ones(n_catchment_eval,1)); % repeat each row of input_thetas n_catchment_eval times
    
    % select random catchment for each model instance
    % iterate over para_sets, each contains n_catchment_eval number of a input_theta
    selected_catchments = zeros(size(input_thetas, 1),1);
    for j = 1:n_para_sets
        ind_s = (j-1)*n_catchment_eval + 1;  % starting index in selected_catchments
        ind_e = j*n_catchment_eval; % end index in selected_catchments
        selected_catchments(ind_s:ind_e,:) = randsample(533,n_catchment_eval); % sample n_catchment_eval for each input_theta
    end
    
    % put the parameters and the selected link into eval_grid for the
    % current model class; size(input_thetas, 1) is the number of instance
    % associated with each class
    starting_ind = (i-1)*size(input_thetas, 1);
    for j=1:size(input_thetas, 1)
        
        ind = j+starting_ind;
        
        eval_grid{ind,1} = model_class; % model class name
        eval_grid{ind,2} = input_thetas(j,:); % model parameters
        eval_grid{ind,3} = selected_catchments(j); % catchment ID
        eval_grid{ind,4} = 0; % evaluation result
    end
end

%% evaluation
tic
step_size = size(input_thetas, 1); % each model class is associated with 'step_size' numbers of evaluation
% iterate over modle class
for i =1:n_model_class
    
    starting_ind = (i-1)*step_size; % index in eval_grid for storing results
    
    % Display progress
    % ---------------------------------------------------------------------
    disp(['model_id = ', num2str(i),' ']);
    disp(['model = ', eval_grid{starting_ind+1,1}]);  % model name
    disp(' ');
    
    % Extract input variables for the current model class
    % ---------------------------------------------------------------------
    model_class = eval_grid{starting_ind+1,1}; % model class
    input_thetas = eval_grid((starting_ind+1):(i*step_size),2); % model parameters
    input_thetas = cell2mat(input_thetas); % cell to matrix
    catchment_ids = eval_grid((starting_ind+1):(i*step_size),3); % catchment_ids
    catchment_ids = cell2mat(catchment_ids);
    
    % parallel model simulation
    % ---------------------------------------------------------------------
    evl_results = zeros(step_size, 3);
    parfor j=1:step_size
        
        % parameters for the model instance
        input_theta = input_thetas(j,:);
        
        % extract catchment data
        catchment_id = catchment_ids(j);
        catchment_data = extract_catchment_data(catchment_datas, catchment_id);
        
        % start model simulation
        [KGE, NSE, RMSE] = catchment_similarity(catchment_data, model_class, input_theta);
        
        % save simulation results to eval_grid
        evl_results(j,:) = [KGE, NSE, RMSE];
    end
    
    % put simulation results to eval_grid
    % ---------------------------------------------------------------------
    for j=1:step_size
        eval_grid{j+starting_ind,4} = evl_results(j,:);
    end
    
    % save the result once evalution for a model clas finishes
    save("sparse_out.mat", "eval_grid")
end
toc

