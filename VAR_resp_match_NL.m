%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code to find Initial Guess for GMA parameters by calculating the GMA model parameters that best fit the VAR-based IRFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

setup.horizon=setup.lags; %horizon up to which IRFS are matched

% for kk=1:setup.horizon+1
%     setup.store_responses(:,:,kk)=VAR_MA_implied_higher_order(setup, kk-1, B_matrix);
% end

%making the code backwards compatible
store_responses=setup.store_responses;

options = optimset('Display','off','TolFun',1e-20,'TolX',1e-20,'MaxIter',500, 'MaxFunEvals',500000);

%first optimize for only 1 Gaussian

for j=1:setup.size_obs
    for k=1:setup.size_obs
        if k>j %above the diagonal
            opt_restr=@(params)var_quad_restr( params,setup,store_responses,j,k ) ;
            %        [xestimate,functionvalue1]=fminsearch(opt_restr,ones(3,1),options);
            [xestimate,functionvalue1]=fminsearch(opt_restr,[.5; 5; 100],options);
            contemp(j,k)=0;
            beta_gen(j,k)=xestimate(1);
            b_gen(j,k)=xestimate(2);
            c_gen(j,k)=xestimate(3);
        elseif j==k
            opt_free=@(params)var_quad_free( params,setup,store_responses,j,k ) ;
            %        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
            [xestimate,functionvalue1]=fminsearch(opt_free,[1 2 1 100],options);
            
            % NOTE: to fine tune the initial guess IG (and also have an IG consistent with some parameter restrictions like b>0 or c>c_lowerbar)
            % we can use fmincon to restrict the IG for a,b and c)
            A=zeros(4,4);
            A(3,3)=-1; %b>-1
            A(4,4)=-1; %b>-1
            b=[0 0 5 -8]';
            
            % fmincon under Ax<=b
            [xestimate,functionvalue1]=fmincon(opt_free,[1 2 1 100],A,b);
            
            beta_diag(j)=(xestimate(1));
            beta_gen(j,k)=xestimate(2);
            b_gen(j,k)=xestimate(3);
            c_gen(j,k)=xestimate(4);
            
        elseif k<j
            opt_below_diag=@(params)var_quad_below_diag( params,setup,store_responses,j,k ) ;
            %        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
            [xestimate,functionvalue1]=fminsearch(opt_below_diag,5*ones(3,1),options);
            
            beta_gen(j,k)=xestimate(1);
            b_gen(j,k)=xestimate(2);
            c_gen(j,k)=xestimate(3);
        end
    end
end

%then optimize for two Gaussians

for j=1:setup.size_obs
    for k=1:setup.size_obs
        if k>j %above the diagonal
            opt_restr=@(params)var_quad_restr_2_gaussian( params,setup,store_responses,j,k ) ;
            %        [xestimate,functionvalue1]=fminsearch(opt_restr,ones(3,1),options);
            [xestimate,functionvalue1]=fminsearch(opt_restr,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 2*b_gen(j,k);c_gen(j,k)],options);
            contemp12(j,k)=0;
            beta_gen12(j,k)=xestimate(1);
            b_gen12(j,k)=xestimate(2);
            c_gen12(j,k)=xestimate(3);
            beta_gen2(j,k)=xestimate(4);
            b_gen2(j,k)=xestimate(5);
            c_gen2(j,k)=xestimate(6);
            
        elseif j==k
            
            opt_free=@(params)var_quad_free_2_gaussian( params,setup,store_responses,j,k ) ;
            %        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
            [xestimate,functionvalue1]=fminsearch(opt_free,[beta_diag(j);beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 2*b_gen(j,k);c_gen(j,k)],options);
            
            beta_diag12(j)=(xestimate(1));
            
            beta_gen12(j,k)=xestimate(2);
            b_gen12(j,k)=xestimate(3);
            c_gen12(j,k)=xestimate(4);
            beta_gen2(j,k)=xestimate(5);
            b_gen2(j,k)=xestimate(6);
            c_gen2(j,k)=xestimate(7);
            
        elseif k<j
            
            opt_below_diag=@(params)var_quad_below_diag_2_gaussian( params,setup,store_responses,j,k ) ;
            %        [xestimate,functionvalue1]=fminsearch(opt_free,ones(4,1),options);
            [xestimate,functionvalue1]=fminsearch(opt_below_diag,[beta_gen(j,k); b_gen(j,k);c_gen(j,k);beta_gen(j,k); 2*b_gen(j,k);c_gen(j,k)],options);
            
            beta_gen12(j,k)=xestimate(1);
            b_gen12(j,k)=xestimate(2);
            c_gen12(j,k)=xestimate(3);
            beta_gen2(j,k)=xestimate(4);
            b_gen2(j,k)=xestimate(5);
            c_gen2(j,k)=xestimate(6);
        end
        
    end
end



%We just take the unrestricted initial response as starting value for the
%unrestricted elements of the contemporaneous matrix
beta_diag=squeeze(setup.store_responses(1:end,setup.index_unrestricted,1));
beta_diag=beta_diag(:);
beta=[];
b=[];
c=[];


%right now only works for up to three Gaussians
beta=[beta;(setup.num_gaussian==1).*beta_gen(:,setup.index_unrestricted)+(setup.num_gaussian==2).*beta_gen12(:,setup.index_unrestricted)];
b=[b;(setup.num_gaussian==1).*b_gen(:,setup.index_unrestricted)+(setup.num_gaussian==2).*b_gen12(:,setup.index_unrestricted)];
c=[c;(setup.num_gaussian==1).*c_gen(:,setup.index_unrestricted)+(setup.num_gaussian==2).*c_gen12(:,setup.index_unrestricted)];
% if max(setup.num_gaussian)>1
beta=[beta;(setup.num_gaussian==2).*beta_gen2(:,setup.index_unrestricted)];
b=[b;(setup.num_gaussian==2).*b_gen2(:,setup.index_unrestricted)];
c=[c;(setup.num_gaussian==2).*c_gen2(:,setup.index_unrestricted)];



beta_diag2=squeeze(setup.store_responses(:,:,1));
beta_diag2=beta_diag2(:);
beta2=[];
b2=[];
c2=[];

setup.num_gaussian2=repmat(setup.num_gaussian,1,setup.size_obs);

%right now, code only works for up to three Gaussians
beta2=[beta2;(setup.num_gaussian2==1).*beta_gen(:,:)+(setup.num_gaussian2==2).*beta_gen12(:,:)];

b2=[b2;(setup.num_gaussian2==1).*b_gen(:,:)+(setup.num_gaussian2==2).*b_gen12(:,:)];
c2=[c2;(setup.num_gaussian2==1).*c_gen(:,:)+(setup.num_gaussian2==2).*c_gen12(:,:)];
% if max(setup.num_gaussian)>1
beta2=[beta2;(setup.num_gaussian2==2).*beta_gen2(:,:)];
b2=[b2;(setup.num_gaussian2==2).*b_gen2(:,:)];
c2=[c2;(setup.num_gaussian2==2).*c_gen2(:,:)];

%beta22=beta;
beta22=beta(beta~=0);
%b22=b;
b22=b(b~=0);
%c22=c;
c22=c(c~=0);
beta2=beta2(beta2~=0);
b2=b2(b2~=0);
c2=c2(c2~=0);
beta_temp=reshape(beta,setup.size_obs,2);
beta_temp=sum(beta_temp,2);

%starting value for optimizers (which in turn return starting value for MCMC)
%setup.initial_parameter=[beta_diag2;beta2(:);b2(:);c2(:);beta_diag;beta22(:);b22(:);c22(:)];
setup.initial_parameter=[beta_diag2;beta2(:);b2(:);c2(:);beta_diag;beta22(:);b22(:);c22(:)];
% Plots the IRFs from a VAR and the fitted GMA impulse responses
check_plots
