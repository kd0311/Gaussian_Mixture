%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=10000;%original is 10000
%number of draws for choosing the scaling matrix in case of standard RW proposal
setup.scaling_draws=2000;%10000
%scaling is recomputed after check_scaling draws
setup.check_scaling=100;%200
%display every disp iter draw
setup.disp_iter=100;
% keep every keep_draw th draw
setup.keep_draw=10;
%proposal=1 ->standard RW
%proposal=2 ->adaptive RW
setup.proposal=1;
%log likelihood computation
%likelihood=3 -> user-supplied logL (only option available in this code)
setup.likelihood=3;
setup.skip_opt=0; %skip optimization and go directly to MCMC

% Original-----> Assign blocks
% setup.number_blocks=4;
% % setup.index_block{1}=[1:2 5 19]';%psi_0
% % setup.index_block{2}=[3 4 6 20]';
% setup.index_block{1}=[1:4 17:18]'; %psi_0
% setup.index_block{2}=[5:8 19:20]'; %a
% setup.index_block{3}=[9:12 21:22]'; %b
% setup.index_block{4}=[13:16 23:24]'; %c


setup.number_blocks=5;
setup.index_block{1}=[3 17]';
setup.index_block{2}=[1 2 4 18]';
setup.index_block{3}=[5:8 19:20]';
setup.index_block{4}=[9:12 21:22]';
setup.index_block{5}=[13:16 23:24]';


%initial scaling for standard RW proposal
%setup.initial_scaling=[.01 .01 .1 .1]'; %one for each block
setup.initial_scaling=[1 200 2 20 20]'; 

%%
% Basic setup for data input
setup.lags=36; %lag length of IRF
setup.size_obs=2; %number of observables
% setup.freq=2;% frequency of data: 1 for monthly, 2 for quarterly
setup.shocks=0; %0-> initial shocks are zero, 1-> initial shocks from VAR (reduces sample size)
setup.polynomials=0; %degree of polynomial detrending
% setup.symmetry=1000; %making sure non of the options are turned on
%setup.VARsym_order=4; %order of symmetric VAR used for starting values
setup.VARsym_order=1; %order of symmetric VAR used for starting values

%setting up symmetry/asymmetry restrictions 
setup.index_restricted=1; %shocks with symmetric IRFs (info)
setup.index_unrestricted=2; %shocks with asymmetric IRFs (mp)

%number of gaussians per IR (same order as observables)
setup.num_gaussian=[1 1]'; %[nb for IRs to shock1, nb for IRs to shock2]
%%
% impose diagonal coefficients of impact matrix to be positive
% beta_1>0 and beta_2>0 and beta_3>0
setup.length_log=3;
setup.index_log=[3 6 20];
setup.length_log=length(setup.index_log);

setup.length_logit_general=0;
setup.index_logit_general=[];
setup.logit_general_lb=[]';
setup.logit_general_ub= []';
%%
setup.length_logit=0;
setup.index_logit=[];

%% Comput IRFs based on LP
% Store LP IRFs
setup.horizon=setup.lags; %horizon up to which IRFS are matched (for initialization and setting of data-driven priors)
setup.store_responses = zeros(2, 2, 36);

%add one more dimension for irf
YSimBar_cbi = reshape(DEU_cbi, [2, 1, 36]);
YSimBar_mp = reshape(DEU_mp, [2, 1, 36]);
setup.store_responses(:, 1, :) = YSimBar_cbi;
setup.store_responses(:, 2, :) = YSimBar_mp;

% for kk=1:setup.horizon
% setup.store_responses(:,:,kk)=VAR_MA_implied_higher_order( setup,kk-1,B_matrix); %store IRFs from VAR
% %setup.store_responses(:, :, kk)=VAR_MA_implied_higher_order(setup, kk-1, B_gdp, B_inflation);
% end

%should additional matrices be stored
setup.add_matrices=1;
setup.sample_size=length(data);
%dimension of those variables
setup.dim_add_matrices=[setup.size_obs setup.sample_size+setup.lags];

%this option is needed because the code would in general allow the
%likelihood function to return estimates of the state if the Kalman filter
%was used
setup.state_size=1;


