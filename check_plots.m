%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code plots the IRFs from a VAR and the fitted GMA impulse responses
% (along with the corresponding Gaussian basis functions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variable_names = {'IP', 'Inflation'};
shock_names = {'CBI Shock','MP Shock'};
%Plot fitted and VAR IRFs
%epsilon_vec=epsilon_vec(1:setup.size_obs,1:setup.lags+1);
epsilon_vec=zeros(setup.size_obs,setup.lags);
 [ Sigma] = unwrap_NL_IRF( setup.initial_parameter,epsilon_vec,setup);
for kk=1:setup.size_obs
    figure;
    for jj=1:setup.size_obs
        subplot(setup.size_obs,1,jj)
        plot(1:size(Sigma,3),squeeze(Sigma(jj,kk,:)),1:size(Sigma,3),squeeze(setup.store_responses(jj,kk,:)),'LineWidth',2)
        legend('fitted GMA response','VAR response')
         grid on
        str = sprintf('%s , %s', variable_names{jj}, shock_names{kk});
        title(str);


    end
end

% Plot the individual Gaussian basis functions
[ Sigma_ind] = unwrap_NL_IRF_individual_bases( setup.initial_parameter,epsilon_vec,setup);

for kk=1:setup.size_obs
    figure;
    for jj=1:setup.size_obs
        subplot(setup.size_obs,1,jj)
        for bb=1:max(setup.num_gaussian) %note that if one irf has less than the maximum number of gaussians, some zero lines will be plotted
            hold on, plot(1:size(Sigma,3),squeeze(Sigma_ind(jj,bb,kk,:)),'LineWidth',2)
        end

        grid on
        str = sprintf('%s , %s', variable_names{jj}, shock_names{kk});
        ylabel(str);
        if jj==1
            title('Individual Gaussian basis functions')
        end

    end
end
%%------------------------Single Chart-----------------------------------------%
% [ Sigma, intercept] = unwrap_NL_IRF( setup.initial_parameter,epsilon_vec,setup);
% for kk = 1:setup.size_obs % Loop over shocks
%     for jj = 1:setup.size_obs % Loop over variables
%         figure; % Create a new figure for each variable-shock combination
%         plot(1:size(Sigma, 3), squeeze(Sigma(jj, kk, :)), ...
%              1:size(Sigma, 3), squeeze(setup.store_responses(jj, kk, :)), 'LineWidth', 2);
%         legend('Fitted GMA response', 'VAR response');
%         grid on;
%         str = sprintf('%s, %s', variable_names{jj}, shock_names{kk});
%         title(str);
%         xlabel('Periods');
%         ylabel('Response');
%     end
% end
% 
% [ Sigma_ind ] = unwrap_NL_IRF_individual_bases( setup.initial_parameter, epsilon_vec, setup);
% 
% output_file = 'Gaussian_Basis_Functions.xlsx';
% 
% % Loop over variables and shocks to generate individual plots and export data
% for kk = 1:setup.size_obs % Loop over shocks
%     for jj = 1:setup.size_obs % Loop over variables
%         figure; % Create a new figure for each variable-shock combination
% 
%         data_to_export = zeros(size(Sigma, 3), max(setup.num_gaussian));
% 
%         for bb = 1:max(setup.num_gaussian) % Loop over Gaussian basis functions
%             hold on;
%             plot(1:size(Sigma, 3), squeeze(Sigma_ind(jj, bb, kk, :)), 'LineWidth', 2);
% 
%             data_to_export(:, bb) = squeeze(Sigma_ind(jj, bb, kk, :));
%         end
% 
%         grid on;
%         str = sprintf('%s, %s', variable_names{jj}, shock_names{kk});
%         ylabel(str);
%         title(sprintf('Individual Gaussian basis functions for %s, %s', variable_names{jj}, shock_names{kk}));
%         xlabel('Periods');
% 
% 
%         sheet_name = sprintf('%s_%s', variable_names{jj}, shock_names{kk});
% 
%         writematrix(data_to_export, output_file, 'Sheet', sheet_name);
% 
%         time_data = (1:size(Sigma, 3))';
%         data_with_time = [time_data, data_to_export];
% 
%         writematrix(data_with_time, output_file, 'Sheet', sheet_name, 'Range', 'A1');
%     end
% end

