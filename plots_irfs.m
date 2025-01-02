%This code plots the estimated responses with error bands

%draws=draws(:,2000:end); %if want to get rid of first N draws as burn-in
setup.number_of_draws=setup.keep_draw*size(draws,2);

Ndraw=100;

inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,Ndraw,1);

% Extract paramters
a_neg_all = zeros(length([9:10]), Ndraw);
b_neg_all = zeros(length([13:14]), Ndraw);
c_neg_all = zeros(length([17:18]), Ndraw);

a_pos_all = zeros(length([21:22]), Ndraw);
b_pos_all = zeros(length([23:24]), Ndraw);
c_pos_all = zeros(length([25:26]), Ndraw)

% Get IRFs to negative shocks

for kk=1:length(indices_for_draws)
    size_of_shock=-1;

    shock_contemp=zeros(setup.size_obs,1);
    shock_contemp(setup.index_unrestricted)=size_of_shock;
    ind_vec=[zeros(1,setup.lags)];
    for jj=1:setup.lags+1
        epsilon=[zeros(setup.size_obs,jj-1) shock_contemp zeros(setup.size_obs,setup.lags+1-jj)];
        [ Sigma, intercept] = unwrap_NL_IRF( draws(:,indices_for_draws(kk)),epsilon,setup);

        IRFs(:,jj,kk)=Sigma(:,setup.index_unrestricted,jj);
    end
    clear ind_vec
    % Extract negative shock parameter values
    a_neg_all(:, kk) = draws([9:10], indices_for_draws(kk));
    b_neg_all(:, kk) = draws([13:14], indices_for_draws(kk));
    c_neg_all(:, kk) = draws([17:18], indices_for_draws(kk));
end

IRFsmed=squeeze(median(IRFs,3));
IRFslow=squeeze(prctile(IRFs,5,3));
IRFshigh=squeeze(prctile(IRFs,95,3));


figure;
for jj=1:setup.size_obs
    subplot(setup.size_obs,1,jj)
    plot(1:setup.lags+1,IRFsmed(jj,:),1:setup.lags+1,IRFslow(jj,:),'k--',1:setup.lags+1,IRFshigh(jj,:),'k--', ...
        1:setup.lags+1,squeeze(setup.store_responses(jj,setup.index_unrestricted,:)),'Linewidth',2)
    title(sprintf('IRF, shock size = %i', size_of_shock))
end
legend('median','5th percentile','95th percentile','VAR')

% print -depsc


% Get IRFs to positive shocks
clear IRFs

inds=setup.number_of_draws/setup.keep_draw;
indices_for_draws=unidrnd(inds,Ndraw,1);


for kk=1:length(indices_for_draws)
    size_of_shock=1;

    shock_contemp=zeros(setup.size_obs,1);
    shock_contemp(setup.index_unrestricted)=size_of_shock;
    ind_vec=[zeros(1,setup.lags)];
    for jj=1:setup.lags+1
        epsilon=[zeros(setup.size_obs,jj-1) shock_contemp zeros(setup.size_obs,setup.lags+1-jj)];
        [ Sigma, intercept] = unwrap_NL_IRF( draws(:,indices_for_draws(kk)),epsilon,setup);

        IRFs(:,jj,kk)=Sigma(:,setup.index_unrestricted,jj);
    end
    clear ind_vec

     % Extract positive shock parameter values
    a_pos_all(:, kk) = draws([21:22], indices_for_draws(kk));
    b_pos_all(:, kk) = draws([23:24], indices_for_draws(kk));
    c_pos_all(:, kk) = draws([25:26], indices_for_draws(kk));
end

IRFsmed=squeeze(median(IRFs,3));
IRFslow=squeeze(prctile(IRFs,5,3));
IRFshigh=squeeze(prctile(IRFs,95,3));


figure;

for jj=1:setup.size_obs
    subplot(setup.size_obs,1,jj)
    plot(1:setup.lags+1,IRFsmed(jj,:),1:setup.lags+1,IRFslow(jj,:),'k--',1:setup.lags+1,IRFshigh(jj,:),'k--',...
        1:setup.lags+1,squeeze(setup.store_responses(jj,setup.index_unrestricted,:)),'Linewidth',2)
    title(sprintf('IRF, shock size = %i', size_of_shock))
end
legend('median','5th percentile','95th percentile','VAR')

% print -depsc

clear IRFs

% %%---------------------------------------------------------------------------
% % Negative shocks IRFs
% Ndraw = 100;
% inds = setup.number_of_draws / setup.keep_draw;
% indices_for_draws = unidrnd(inds, Ndraw, 1);
% 
% a_neg_all = zeros(length([9:10]), Ndraw);
% b_neg_all = zeros(length([13:14]), Ndraw);
% c_neg_all = zeros(length([17:18]), Ndraw);
% 
% % Get IRFs to negative shocks
% IRFs = zeros(setup.size_obs, setup.lags + 1, Ndraw);
% for kk = 1:length(indices_for_draws)
%     size_of_shock = -1;
% 
%     shock_contemp = zeros(setup.size_obs, 1);
%     shock_contemp(setup.index_unrestricted) = size_of_shock;
% 
%     for jj = 1:setup.lags + 1
%         epsilon = [zeros(setup.size_obs, jj - 1), shock_contemp, zeros(setup.size_obs, setup.lags + 1 - jj)];
%         [Sigma, intercept] = unwrap_NL_IRF(draws(:, indices_for_draws(kk)), epsilon, setup);
%         IRFs(:, jj, kk) = Sigma(:, setup.index_unrestricted, jj);
%     end
% 
%     % Extract negative shock parameter values
%     a_neg_all(:, kk) = draws([9:10], indices_for_draws(kk));
%     b_neg_all(:, kk) = draws([13:14], indices_for_draws(kk));
%     c_neg_all(:, kk) = draws([17:18], indices_for_draws(kk));
% end
% 
% IRFsmed = squeeze(median(IRFs, 3));
% IRFslow = squeeze(prctile(IRFs, 5, 3));
% IRFshigh = squeeze(prctile(IRFs, 95, 3));
% 
% % Separate figures for each variable (negative shocks)
% for jj = 1:setup.size_obs
%     figure;
%     plot(1:setup.lags + 1, IRFsmed(jj,:), 1:setup.lags + 1, IRFslow(jj,:), 'k--', 1:setup.lags + 1, IRFshigh(jj,:), 'k--', ...
%         1:setup.lags + 1, squeeze(setup.store_responses(jj, setup.index_unrestricted, :)), 'Linewidth', 2);
%     legend('median', '5th percentile', '95th percentile', 'VAR');
% end
% 
% % Positive shocks IRFs
% 
% a_pos_all = zeros(length([21:22]), Ndraw);
% b_pos_all = zeros(length([23:24]), Ndraw);
% c_pos_all = zeros(length([25:26]), Ndraw);
% 
% IRFs = zeros(setup.size_obs, setup.lags + 1, Ndraw); % Clear previous IRFs
% for kk = 1:length(indices_for_draws)
%     size_of_shock = 1;
% 
%     shock_contemp = zeros(setup.size_obs, 1);
%     shock_contemp(setup.index_unrestricted) = size_of_shock;
% 
%     for jj = 1:setup.lags + 1
%         epsilon = [zeros(setup.size_obs, jj - 1), shock_contemp, zeros(setup.size_obs, setup.lags + 1 - jj)];
%         [Sigma, intercept] = unwrap_NL_IRF(draws(:, indices_for_draws(kk)), epsilon, setup);
%         IRFs(:, jj, kk) = Sigma(:, setup.index_unrestricted, jj);
%     end
% 
%     % Extract positive shock parameter values
%     a_pos_all(:, kk) = draws([21:22], indices_for_draws(kk));
%     b_pos_all(:, kk) = draws([23:24], indices_for_draws(kk));
%     c_pos_all(:, kk) = draws([25:26], indices_for_draws(kk));
% end
% 
% IRFsmed = squeeze(median(IRFs, 3));
% IRFslow = squeeze(prctile(IRFs, 5, 3));
% IRFshigh = squeeze(prctile(IRFs, 95, 3));
% 
% % Separate figures for each variable (positive shocks)
% for jj = 1:setup.size_obs
%     figure;
%     plot(1:setup.lags + 1, IRFsmed(jj,:), 1:setup.lags + 1, IRFslow(jj,:), 'k--', 1:setup.lags + 1, IRFshigh(jj,:), 'k--', ...
%         1:setup.lags + 1, squeeze(setup.store_responses(jj, setup.index_unrestricted, :)), 'Linewidth', 2);
%     legend('median', '5th percentile', '95th percentile', 'VAR');
% end
