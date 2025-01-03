function [ VAR_match ] = var_quad_restr( params,setup,store_responses,j,k )
%objective function for VAR responses matching
[ MA_matrices ] = wrap_BM_var_opt_12( params, setup );
MA_h=MA_matrices(1,1,1:setup.horizon);
sr=store_responses(j,k,:);
% size(MA_h)
% size(sr)

VAR_match=(MA_h(:)-sr(:))'*(MA_h(:)-sr(:));
end

