%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main prg to run the FAIR asymmetric estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

working_folder = 'C:\Users\kdai\Desktop\FAIR\FAIR_2Var';
addpath(working_folder);
cd(working_folder);

filedata = [working_folder '\EA20Data.xlsx'];
sheetdata = 'DataM1';
rangedata = 'A1:KM5'; 
[data,~,variables] = xlsread(filedata,sheetdata,rangedata);
%variables(2:end,:) = []; variables(:,1) = [];

%Import variables and 2 external shocks
variables_data = data(1:2, :);  
ZZZ = data(3:4, :);
cbi_shock = data(3, :);  
mp_shock = data(4, :);   
data = variables_data;

% gdp_data = data(1, :); 
% gdp_data = gdp_data'; 
% [h, pValue] = adftest(gdp_data);  
% disp(['ADF test: h = ', num2str(h), ', p-value = ', num2str(pValue)]);

clear pathdata filedata sheetdata rangedata
save 'C:\Users\kdai\Desktop\FAIR\FAIR_2Var\data\data_file2.mat';

tic

mkdir results

addpath(['lib']);
addpath(['results']);
%%
%Initialization/Parametrization of routine
setup_NL;
%%
%To find starting values (or, in our case, priors), match IRFs with those of VAR:
VAR_resp_match_NL;
setup.length_param_vector=length(setup.initial_parameter);

%close all
pause(.1)


%% Set parameter restrictions

% Impose lower-bnd / upper-bnd on parameters
% logit - parameter has to be between 0 and 1
% log - parameter has to be non-negative
% logit_general - upper and lower bounds are set by hand

% % Impose -50<b<50 and -50<c<500 (for illustration)
% setup.index_logit_general=[11:14 27:30 15:18 31:34];
% setup.length_logit_general=length(setup.index_logit_general);
% setup.logit_general_lb=[-50*ones(1,8) -50*ones(1,8)]';
% setup.logit_general_ub=[50*ones(1,8) 500*ones(1,8)]';

% setup.index_logit_general=[11:14 23:24 15:18 25:26];
% setup.length_logit_general=length(setup.index_logit_general);
% setup.logit_general_lb=[-50*ones(1,6) -50*ones(1,6)]';
% setup.logit_general_ub=[50*ones(1,6) 500*ones(1,6)]'

% % Alternative identification scheme: sign restriction to (set) identify shock
% % beta0_pi<0 and a_pi<0 for shock to third variable (monetary shock)
% % a_r>0 (and beta0_r>0, which is always imposed in setup_NL) 
% % also impose -25<b<25 and 8<c<4000 (for illustration)
% setup.index_logit_general=[4 6 8 10 20 22 24 26 11:14 27:30 15:18 31:34];
% setup.length_logit_general=length(setup.index_logit_general);
% setup.logit_general_lb=[-10 -10 -10 -10 -10 -10 -10 -10 -100*ones(1,8) -100*ones(1,8)]';
% setup.logit_general_ub=[10 10 10 10 10 10 10 10 500*ones(1,8) 5000*ones(1,8)]';

 % Impose -100<b<1000 and -50<c<500 (for illustration)
setup.index_logit_general=[6 20 10 22 11:14 23:24 15:18 25:26];
setup.length_logit_general=length(setup.index_logit_general);
setup.logit_general_lb=[-20 -20 -20 -20 -100*ones(1,6) -100*ones(1,6)]';
setup.logit_general_ub=[200 200 200 200 1000*ones(1,6) 1000*ones(1,6)]';

%Verify that initial guess is inside the bounds:
if max(setup.initial_parameter(setup.index_logit_general)<=setup.logit_general_lb)==1 | max(setup.initial_parameter(setup.index_logit_general)>=setup.logit_general_ub)==1
    'Initial Guess inside the bounds!!'
    stop
end

%% Set prior parameters (and identifying restrictions imposed through priors)
% %scaling for adaptive MCMC (see handbook of MCMC, page 104) ADAPTED
setup.scaling_adaptive=.1^2/setup.length_param_vector;

setup.index_gamma=[];
setup.index_normal=[];

%Only use Normal priors, but we could set for gamma priors so the prior
%distribution will be normal+gamma (in prior sub-function)
setup.index_normal=[1:setup.length_param_vector]';
setup.normal_prior_means=[setup.initial_parameter;];
setup.normal_prior_std = 3 * abs([setup.initial_parameter;]);


%contemporaneous matrix:
setup.normal_prior_std(setup.index_block{1})=10;
% a:
setup.normal_prior_std(setup.index_block{2})=10;
% b:
setup.normal_prior_std(setup.index_block{3})=10;
% c:
setup.normal_prior_std(setup.index_block{4})=20^2;

% Example of short-run pointwise restriction: tight prior at 0 for upper diagonal elements of last column of contemp impact matrix)
%setup.normal_prior_std([5 6 19 20])=[1e-3,1e1,1e-3,1e-3];


%scaling for adaptive MCMC (see handbook of MCMC, page 104) ADAPTED
setup.scaling_adaptive=.1^2/setup.length_param_vector;

%% Test parameter
data_mean = mean(data, 2);
data_std = std(data, 0, 2);
disp('Data Mean:');
disp(data_mean);
disp('Data Standard Deviation:');
disp(data_std);

figure;

% hist
subplot(4, 1, 1);
histogram(data(:));
title('Data Distribution');
xlabel('Data Value');
ylabel('Frequency');

% box
subplot(4, 1, 2);
boxplot(data(:));
title('Data Boxplot');
ylabel('Data Value');

% initial param hist
subplot(4, 1, 3);
initial_params = setup.initial_parameter;
histogram(initial_params);
title('Initial Parameter Estimates');
xlabel('Parameter Value');
ylabel('Frequency');
ylim([0, 30]);

% initial param box
subplot(4, 1, 4);
boxplot(initial_params);
title('Initial Parameter Boxplot');
xlabel('Parameter');
ylabel('Value');
%% Estimation
[ draws, acc_rate, log_posteriors, statedraws, add_matrices] = sampling_MH( setup );
%add_matrices store the estimated shocks

save results/test_NL
toc
%% EXPORT a b c
a_neg = draws(9:10, :);
b_neg = draws(13:14, :);
c_neg= draws(17:18, :);
a_pos = draws(21:22, :);
b_pos = draws(23:24, :);
c_pos = draws(25:26, :);

% disp(a_neg);
% disp(b_neg);
% disp(c_neg);
% disp(a_pos);
% disp(b_pos);
% disp(c_pos);

%% Get shock series
inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,1000,1);
shocks=add_matrices(:,:,indices_for_draws);
median_shock=squeeze(prctile(shocks,50,3));
lower_shock=squeeze(prctile(shocks,5,3));
upper_shock=squeeze(prctile(shocks,95,3));

% Plot IRFs
plots_irfs


%%%%%%%%% Alternative MLE 
%% Likelihood Function for MLE
%function log_l = neg_likelihood(params, data, setup)
    %[log_l, ~, ~, ~] = likelihood(data, params, setup);
    %log_l = -log_l; % Negate cuz we will minimize
%end

%% Perform MLE Estimation
%options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxIter', 1000, 'Display', 'iter', 'MaxFunEvals', 10000);
%[param_estimates, fval, exitflag, output] = fminunc(@(params) neg_likelihood(params, data, setup), setup.initial_parameter, options);

%save('results/estimation_MLE.mat', 'param_estimates', 'fval', 'exitflag', 'output');
%toc

% Get estimated shocks
%setup.initial_parameter = param_estimates;
%[~, ~, epsilon, uvec] = likelihood(data, param_estimates, setup);

%%
% Get shock series
%shocks = epsilon;
%median_shock = median(shocks, 2);
%lower_shock = prctile(shocks, 5, 2);
%upper_shock = prctile(shocks, 95, 2);

% Plot IRFs
%plots_irfs;
