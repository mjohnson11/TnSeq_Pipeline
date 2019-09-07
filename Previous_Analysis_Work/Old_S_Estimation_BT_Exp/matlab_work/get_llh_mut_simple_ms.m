function llh = get_llh_mut_simple_ms( X, mu, ptr, sigerr, mutr, sigtr)
% GET_LLH_MUT_SIMPLE_MS calulate the negative of the log-likelihood of data
% under the model with outliers
%
% LLH = GET_LLH_MUT_SIMPLE_MS(X, MU, PTR, SIGERR, MUTR, SIGTR) calculates
% the likelihood of observing a vector of barcode fitnesses X (one
% measurement per barcode), given the parameters MU (true selection
% coefficient of mutation), PTR (probability of a transformation artifact),
% SIGERR (variance of the error distribution), MUTR (mean effect of a
% transformation artifact), SIGTR (variance of the transformation
% artifact). For the details of the model see Johnson, Martsul,
% Kryazhimskiy, Desai (2019). Higher fitness yeast genotypes are less
% robust to deleterious mutation. LLH is the value of the negative of the
% log-likelihood

logthr = 100; % threshold on the difference between the logs of probabilities 

n_bc = size(X,1);

llh = 0;

for ibc = 1:n_bc
        
    % This is the term that corresponds to the absence of a
    % transformation artifact
    A = logmvn( X(ibc) ,  mu , sigerr );
    
    % This is the term that corresponds to the presence of a
    % transformation artifact
    if ptr > 0
        B = logmvn( X(ibc) ,  mu - mutr, sigerr + sigtr);
    end

    % Summing the two terms, noting that we already actually calculated
    % logs of the probabilities:
    if ptr == 0
        llh = llh + A;
    elseif ptr == 1
        llh = llh + B;
    else
        
        
        if A - B < -logthr
        % if A term is much smaller than the B term, only retain B            
            llh = llh + B + log(ptr);
        
        elseif A - B > logthr
        % if B term is much smaller than the A term, only retain A
            llh = llh + A + log(1-ptr);
            
        else
        % otherwise retain both
            llh = llh + log( (1-ptr)*exp(A) + ptr*exp(B) );
        end
    end
end
llh = -llh;



    % log of the normal probability density function
    function r = logmvn(x, m, sig)
        d = size(m,1);
        r = -1/2 * (d * log(2*pi) + log(det(sig)) ) - 1/2 * (x - m)' * (sig \ (x-m));
    end


end