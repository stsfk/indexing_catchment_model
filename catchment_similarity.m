function [KGE, NSE, RMSE] = catchment_similarity(catchment_data,model,input_theta)
%catchment_similarity

%% Define the solver settings
% Create a solver settings data input structure.
% NOTE: the names of all structure fields are hard-coded in each model
% file. These should not be changed.
input_solver.name              = 'createOdeApprox_IE';        % Use Implicit Euler to approximate ODE's
input_solver.resnorm_tolerance = 0.1;                         % Root-finding convergence tolerance
input_solver.resnorm_maxiter   = 6;

%% data precessing
% warm-up data
warm_up_catchment_data = catchment_data(1:366,:);             % 1-year data

warm_up_climatology.precip   = warm_up_catchment_data.P;      % Daily data: P rate  [mm/d]
warm_up_climatology.temp     = warm_up_catchment_data.T;      % Daily data: mean T  [degree C]
warm_up_climatology.pet      = warm_up_catchment_data.PET;    % Daily data: Ep rate [mm/d]
warm_up_climatology.Q        = warm_up_catchment_data.Q;      % Daily data: Q rate [mm/d]
warm_up_climatology.delta_t  = 1;

% main forcing data
main_data = catchment_data(367:4018,:);                        % From the second year to the end

main_climatology.precip   = main_data.P;                      % Daily data: P rate  [mm/d]
main_climatology.temp     = main_data.T;                      % Daily data: mean T  [degree C]
main_climatology.pet      = main_data.PET;                    % Daily data: Ep rate [mm/d]
main_climatology.Q        = main_data.Q;                      % Daily data: Q rate [mm/d]
main_climatology.delta_t  = 1;

%% warm-up
numStore    = str2double(model(end-1));
input_s0    = zeros(numStore,1);

% First iteration, storage starts from 0
[~,...                                            % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
    ~,...                                         % Internal model fluxes
    output_ss] = ...                              % Internal storages
    feval(model,...                               % Model function name
    warm_up_climatology,...                       % Time series of climatic fluxes in simulation period
    input_s0,...                                  % Initial storages
    input_theta,...                               % Parameter values
    input_solver);

input_s0 = structfun(@last_element,output_ss);        % Storage level prior to simulation

thd = 0.05;                                           % Threshold of storage change ratio
for iter=2:20
    
    [~,...                                            % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
        ~,...                                         % Internal model fluxes
        output_ss] = ...                              % Internal storages
        feval(model,...                               % Model function name
        warm_up_climatology,...                       % Time series of climatic fluxes in simulation period
        input_s0,...                                  % Initial storages
        input_theta,...                               % Parameter values
        input_solver);
    
    input_s0_new = structfun(@last_element,output_ss);                  % new storage to be use in next iteration
    storage_change = abs(input_s0_new - input_s0);                      % storage change absolute value
    
    if max(storage_change - thd * abs(input_s0)) <= 0                   % compare to thd storage change threshold, abs to account for deficit store
        break
    end
    
    input_s0 = input_s0_new;                                            % use storage at the end as initial storage for next step
end

%% Main period
% Set the inital storages
input_s0  = structfun(@last_element,output_ss);

[output_ex,...                                    % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
    ~,...                                         % Internal model fluxes
    ~] = ...                                      % Internal storages
    feval(model,...                               % Model function name
    main_climatology,...                          % Time series of climatic fluxes in simulation period
    input_s0,...                                  % Initial storages
    input_theta,...                               % Parameter values
    input_solver);

%% Evaluation
KGE = of_KGE(main_climatology.Q, output_ex.Q);
NSE = of_NSE(main_climatology.Q, output_ex.Q);
RMSE = of_RMSE(main_climatology.Q, output_ex.Q);

end

