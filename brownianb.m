function pdf = brownianb(d, step, tol, max_it)

    if (nargin < 4)
        max_it = 1e3;  % Maximum iterations for sigma_m steepest descent 
    end
    if (nargin < 3)
        tol = 1e-3;    % Threshold for gradient magnitude in sigma_m s.d.
    end
    if (nargin < 2)
        step = 0.2;    % Learning step for sigma_m s.d.
    end
    if (nargin < 1)
        % This is a recorded trajectory example
        d.x = [1,2; 3,4; 5,5; 7,8; 8,7; 7,6; 6,5]; % coordinates
        d.t = [0,   1,   2,   3,   4,   5,   6  ]; % observation times
        d.s = 0.1;        % observation error (can be global or per-sample)
    end
    
    % The algorithm always considers distinct measure errors per sample
    if (isscalar(d.s))
        d.s = d.s * ones(size(d.t));
    end

    % We have n observations
    n_obs = length(d.t);

    % We construct bridges skipping data: 1=3, 3=5, ...
    % and leave the intermediates out for fitting sigma_m with ML
    
    % List of indices of data points left for fitting ML
    % Note we cannot fit the last one because there is no bridge
    idx_fitting = (2:2:n_obs-1)';
    
    % List of begin and end indices for each training bridge
    idx_model = 1:2:n_obs;
    bridges = [idx_model(1:end-1); idx_model(2:end)]';
    cached_bridge_t_end = d.t(bridges(:,2));

    % We have to fit this parameter
    sigma_m = 1; % start with unit variance
    
    % Finds mean, variance and alpha of expected position at time t
    function [mu, sigma, alpha] = stats(t,i)
        % find current bridge
        if (nargin < 2)
            i = search(t, cached_bridge_t_end);
            i_begin = bridges(i,1);   % begin index of the current bridge
            i_end = bridges(i,2);     % end index of the current bridge
        else
            i_begin = i(:, 1);
            i_end = i(:, 2);
        end
        dur = d.t(i_end) - d.t(i_begin); % duration of the current bridge 
        pos = (t - d.t(i_begin)) ./ dur;  % time position within the bridge
        % mean of expected particle position at time t
        pp = repmat(pos, [1 2]);
        mu = (1-pp) .* d.x(i_begin,:) + pp .* d.x(i_end,:);
        % alpha is dur*pos*(1-pos), for convenience
        % (it's the derivative of sigma wrt sigma_m)
        alpha = dur.*pos.*(1-pos);
        % variance of expected particle position at time t 
        sigma = alpha*sigma_m + (1-pos).^2.*d.s(i_begin) + pos.^2.*d.s(i_end);
    end

    % Derivative of the negative log likelihood of the model (given sigma_m)
    % which is the sum of log(p_xt) for each (x,t) in the fitting group
    % see paper for details
    % negative log-likelihood without bias and offset added for debugging
    function [dnll,nll] = dnloglike()
        dnll = 0;
        nll = 0;
        for i = 1:length(idx_fitting)
            idx = idx_fitting(i);
            % statistics of where the current fitting point should be
            % according to the current model
            [mu,sigma,alpha] = stats(d.t(idx), [idx-1,idx+1]);
            % derivative of the neg log-likelihood of current fitting point
            % is aggregated to the total
            dist = d.x(idx,:) - mu;
            dd = dist*dist'/(2*sigma);
            dnll = dnll + alpha / sigma * (1 - dd);
            % negative log-likelihood, for debugging
            % nll = nll + log(sigma) + dd;
        end
    end

    % Optimal sigma_m is found by minimizing negative-log-likelihood
    % Using a simple and naive steepest descend
    % we do not know if the function is convex or not!
    dnll = Inf;
    it = 0;
    while (abs(dnll) > tol && it < max_it)
        [dnll,~] = dnloglike();
        % disp([nll, dnll, sigma_m]);
        sigma_m = max(sigma_m - step * dnll, eps);
        it = it + 1;
    end
    
    % The result of this algorithm is the numerical pdf of the particle
    % position x at a given instant t
    sq2p = 1/sqrt(2*pi);
    function p = p_xt(x,t)
        sz = size(t);
        t = t(:);
        % lookup on the normal distribution
        [m,s] = stats(t);
        xm = bsxfun(@minus, x, m);
        p = sq2p./s .* exp(-sum(xm.*xm,2)./(s*2));
        p = reshape(p, sz);
    end
    pdf = @p_xt;

end