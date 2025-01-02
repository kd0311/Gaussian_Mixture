%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%General Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%number of draws after burn-in 
setup.number_of_draws=50000;%original is 10000
%number of draws for choosing the scaling matrix in case of standard RW proposal
setup.scaling_draws=4000;%10000
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

% Assign blocks
setup.number_blocks=4;
% setup.index_block{1}=[1:2 5 19]';%psi_0
% setup.index_block{2}=[3 4 6 20]';
setup.index_block{1}=[1:6 19:20]'; %psi_0
setup.index_block{3}=[7:10 21:22]'; %a
setup.index_block{4}=[11:14 23:24]'; %b
setup.index_block{5}=[15:18 25:26]'; %c


% setup.number_blocks=1;
% setup.index_block{1}=[1:51]';
%initial scaling for standard RW proposal
setup.initial_scaling=[.01 .01 .1 .1]'; %one for each block
%setup.initial_scaling=[1 200 2 20 20]'; 

%%
% Basic setup for data input
setup.lags=36; %lag length of IRF
setup.size_obs=2; %number of observables
% setup.freq=2;% frequency of data: 1 for monthly, 2 for quarterly
setup.shocks=0; %0-> initial shocks are zero, 1-> initial shocks from VAR (reduces sample size)
setup.polynomials=0; %degree of polynomial detrending
% setup.symmetry=1000; %making sure non of the options are turned on
%setup.VARsym_order=1; %order of symmetric VAR used for starting values
setup.VARsym_order=4; %order of symmetric VAR used for starting values

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

% Load the data file, add polynomial trend if needed and estimate VAR
load (['data/data_file',num2str(setup.size_obs),'.mat'])
run (['data/detrending']); % include polynomial detrending if needed
%%
setup.data = data;
ip = data(1,:);
inf = data(2,:);
ip_data = ip';
inf_data = inf';
oosT = 36;
estT = 100 - oosT;
estIdx = 1:estT;     % Estimation sample indices
oosIdx = (100 - 35):100; % Out-of-sample indices
EstY = [ip_data(estIdx) inf_data(estIdx)]; % In-sample responses
estX = cbi_shock(estIdx);       
estX2 = mp_shock(estIdx);   % In-sample exogenous data
nn = size(EstY,2);

OOSY = [ip_data(oosIdx) inf_data(oosIdx)]; % Out-of-sample responses
oosX = cbi_shock(oosIdx); 
oosX2 = mp_shock(oosIdx); 

Mdl = varm(nn,4);
EstMdl = estimate(Mdl,EstY,'X',estX);
EstMdl2 = estimate(Mdl,EstY,'X',estX2);
summarize(EstMdl)
summarize(EstMdl2)
numPaths = 1000;
Y0 = EstY((end-3):end,:);
rng(1); % For reproducibility
YSim_cbi = simulate(EstMdl,oosT,'X',oosX,'Y0',Y0,'NumPaths',numPaths);
YSimBar_cbi = mean(YSim_cbi,3);   % 3 refers to 3rd dimension, which is number of total simulations
YSim_mp = simulate(EstMdl2,oosT,'X',oosX2,'Y0',Y0,'NumPaths',numPaths);
YSimBar_mp= mean(YSim_mp,3);

YSimBar_cbi = YSimBar_cbi';
YSimBar_mp = YSimBar_mp';

