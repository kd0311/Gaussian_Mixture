%% Comput IRFs based on VAR
% Store VAR IRFs
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
