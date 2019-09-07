function myfit = fit_data_norm_ms( X, init, RelTol, kmax )
% FIT_DATA_NORM_MS fits selection coefficients of all mutations in multiple
% strains using maximum likelihood in a model that takes outlier barcodes
% into account
%
% MYFIT = FIT_DATA_NORM_MS( X, INIT, RELTOL, ABSTOL, KMAX ) takes an M by K
% cell array X of the barcode selection coefficients and estimates the
% selection coefficients of mutations, where M is the number of mutations
% and K is the number of strains. X{imut,istrain} is a vector of selection
% coefficients of barcodes that correspond to mutation imut in strain
% istrain. All other input arguments are optional. INIT is the initial
% values for all parameters.
%
% INIT must be a structure with the following
% fields: mumat is an M by K matrix of selection coefficients of
% mutations; ptr is the probability of a transformation artifact; mutr is
% the mean effect of the transformation artifact; sigtr is the variance of
% the transformation artifact; sigerr is the error variance. If these
% values are not specified the following values are used: ptr = 0.05;
% sigerr = 0.05; mutr = 0.2; sigtr = 0.25;
%
% RELTOL is the relative tolerance used for stopping the maximization
% procedure. Default value is 0.001. 
%
% KMAX is the maximum number of iterations to maximize likelihood
%
% MYFIT is a structure with the following fields: mumat is an M by K matrix
% of estimated selection coefficients of mutations; ptr is the estimated
% probability of a transformation artifact; mutr is the estimated mean
% effect of the transformation artifact; sigtr is the estimated variance of
% the transformation artifact; sigerr is the estimated error variance.

[n_mut , n_strain] = size( X );
% n_mut = number of mutations
% n_strain = number of strains

% setting initial values:
if nargin < 2
    curr_mumat = cellfun(@mean, X, 'UniformOutput', true);
    % initialize the selection coefficients of mutations using the means
    % across barcodes
    
    curr_ptr = 0.05;
    curr_sigerr = 0.25;
    curr_mutr = 0.20;
    curr_sigtr = 4;
    
    RelTol = 1e-5;
    kmax = 100;
else
    curr_mumat = init.mumat;
    
    curr_ptr = init.ptr;
    curr_sigtr = init.sigtr;
    curr_mutr = init.mutr;
    curr_sigerr = init.sigerr;

    if nargin < 4
        kmax = 100;
    end
end

AbsToll = 1e-6;

Ntr = cellfun(@length, X); % number of transformants per strain per mutation
NtrTot = nansum( Ntr(:) ); % total number of transformants
NtrMax = nanmax( Ntr(:) ); % maximum number of transformats per strain per mutation

% specifying lower bounds for global parameters:
param_lb.ptr = 0;
param_lb.sigerr = 0;
param_lb.mutr = -Inf;
param_lb.sigtr = 0;

lb = var2vec( param_lb.ptr, param_lb.sigerr, param_lb.mutr, param_lb.sigtr );

% specifying upper bounds for global parameters
param_ub.ptr = 1;
param_ub.sigerr = Inf;
param_ub.mutr = Inf;
param_ub.sigtr = Inf;

ub = var2vec(param_ub.ptr,  param_ub.sigerr,  param_ub.mutr, param_ub.sigtr); 


options = optimset('Display', 'notify');


old_llh = Inf;
IFSTOP = false;

k = 1;
while ~IFSTOP
    
    % First step: optimize the selection coefficients
    for istrain = 1:n_strain
        for imut = 1:n_mut
            
            if isempty( X{imut,istrain} )
                continue;
            end
            
            % maximize the llh starting from mu = 0:
            [x0, llh0] = fminsearch(@(mu) get_llh_mut_simple_ms( X{imut,istrain},...
                mu, curr_ptr, curr_sigerr, curr_mutr, curr_sigtr), 0 );
            
            % maximize the llh starting from current value of the selection coefficient for this mutation:
            [x1, llh1] = fminsearch(@(mu) get_llh_mut_simple_ms(X{imut,istrain}, ...
                mu, curr_ptr, curr_sigerr, curr_mutr, curr_sigtr), curr_mumat(imut,istrain) );

            % If arrived at different points, choose the one that gives higher llh
            if llh0 < llh1
                curr_mumat(imut,istrain,:) = x0;
            else
                curr_mumat(imut,istrain,:) = x1;
            end
        end
    end
    clear x0 llh0 x1 llh1;
    
    
    % Second step: optimize the global parameters ptr, sigerr, mutr, sigtr
    p0 = var2vec(curr_ptr, curr_sigerr, curr_mutr , curr_sigtr);
    [p,curr_llh,exitflag,output] = fmincon(@get_llh_global, p0, [], [], [], [], lb, ub, @mycon, options);    
    [curr_ptr, curr_sigerr, curr_mutr, curr_sigtr] = vec2var( p );

    % Decide whether to stop or not:
    LLH_REL = abs(curr_llh - old_llh)/abs(old_llh) < RelTol;
    PARAM_REL = abs(p(p0>AbsToll) - p0(p0>AbsToll))./p0(p0>AbsToll) < RelTol;    
        
    if  (  LLH_REL && all( PARAM_REL) ) || k >= kmax
        IFSTOP = true;
    end
    
    % Taking special care of the Ptr parameter because it cannot be
    % too small or too large:
    
    if curr_ptr < 1/NtrMax
    % We really should not be able to infer a transformation artifact
    % probability that is lower than the inverse of the maximum number of
    % transformants per strain per mutation
    
        curr_ptr = 0;
        curr_mutr = 0;
        curr_sigtr = 0;
        p1 = var2vec(curr_ptr, curr_sigerr, curr_mutr , curr_sigtr);
        curr_llh = get_llh_global(p1);
        fprintf('ptr = 0: curr_llh = %.6f, adj_llh = %.6f\n', get_llh_global(p), curr_llh);

    elseif 1-curr_ptr  < 1/NtrMax
    % On the other hand, if ptr approaches 1, the model simplifies to a
    % model without transformation artifacts but with larger error
    % variance and all selection coefficients shifted by mutr
    
        curr_ptr = 1;
        p1 = var2vec(curr_ptr, curr_sigerr, curr_mutr , curr_sigtr);
        curr_llh = get_llh_global(p1);
        fprintf('ptr = 1: curr_llh = %.6f, adj_llh = %.6f\n', get_llh_global(p), curr_llh);

        curr_ptr = 0;
        curr_mumat = curr_mumat - curr_mutr;
        curr_mutr = 0;
        curr_sigerr = curr_sigtr + curr_sigerr;
        curr_sigtr = 0;
                
        p1 = var2vec(curr_ptr, curr_sigerr, curr_mutr , curr_sigtr);
        curr_llh = get_llh_global(p1);
        fprintf('=> ptr = 0: curr_llh = %.6f, adj_llh = %.6f\n', get_llh_global(p), curr_llh);
    end
    
    old_llh = curr_llh;   
    fprintf('k = %d, llh = %.6f\n', k, curr_llh);        
                
    k = k+1;
end



%% Third, re-optimize selection coefficients and calculate P-values for selection coefficient

myfit.pvals = nan(size(curr_mumat));

nstart = 5; % number of starting points to try for each mutation in each strain
startmuvec = [0 ,  linspace(nanmin(curr_mumat(:)), nanmax(curr_mumat(:)), nstart-1)] ;

for istrain = 1:n_strain
    for imut = 1:n_mut
        
        if isempty( X{imut,istrain} )
            continue;
        end

        % Here we really attempt to find the right value of the selection
        % coefficient that maximizes llh, so we are starting from multiple
        % starting points
        
        mu_min = startmuvec(:,1);
        llh_min = Inf;
        
        for istart = 1:length(startmuvec)
            [x_curr, llh_curr] = fminsearch(@(mu)...
                get_llh_mut_simple_ms(X{imut,istrain}, mu, curr_ptr, curr_sigerr, curr_mutr, curr_sigtr),...
                startmuvec(istart) );
            
            if llh_curr < llh_min
                llh_min = llh_curr;
                mu_min = x_curr;
            end
        end
        curr_mumat(imut,istrain) = mu_min;
        
        % null model:
        llh0 = get_llh_mut_simple_ms(X{imut,istrain}, 0, curr_ptr, curr_sigerr, curr_mutr, curr_sigtr);
        
        x = -2 * (llh_min - llh0);
        
        if x < 0
            % This should not happen: likelihood under the null model
            % should always be lower
            fprintf('Mut #%d, Env %s: LLH0 = %.6f, LLH1 = %.6f\n',...
                imut, ienv, llh0, llh1);
        end
        myfit.pvals(imut, istrain) = 1 - chi2cdf(x, 1);
    end
end
clear x0 llh0 x1 llh1;


myfit.mutr = curr_mutr;
myfit.sigtr = curr_sigtr;
myfit.ptr = curr_ptr;
myfit.sigerr = curr_sigerr;
myfit.mumat = curr_mumat;









    function llh = get_llh_global( p )
        [ptr, sigerr, mutr, sigtr] = vec2var( p );       
        
        llh = 0;
        for istrain1 = 1:n_strain
            for imut1 = 1:n_mut
                if isempty( X{imut1,istrain1} )
                    continue;
                end
                llh = llh + get_llh_mut_simple_ms( X{imut1,istrain1}, curr_mumat(imut1,istrain1), ptr, sigerr, mutr, sigtr);
            end
        end
    end


    function r = var2vec( ptr,  sigerr,  mutr, sigtr )
        r = [ptr; sigerr; mutr; sigtr];
    end


    function [ptr, sigerr, mutr, sigtr] = vec2var( x )
        ptr = x(1);
        sigerr = x(2);
        mutr = x(3);
        sigtr = x(4);
    end

    function [c, ceq] = mycon(x)
        [ptr, sigerr, mutr, sigtr] = vec2var(x);
        c = -min(eig( sigtr ));
        ceq = 0;
    end
end