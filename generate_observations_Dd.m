function [X, group_elements, rho_emp] = generate_observations_Dd(x, M, sigma, rho, noisetype)
% 
% Given a signal x of length d, generates a matrix X of size d x M such
% that each column of X is a randomly transformed version of x under the
% dihedral group D_d (shifts and reflections), with
% iid Gaussian noise of variance sigma^2 added on top. If x is complex,
% the noise is also complex.

    x = x(:);
    d = length(x);
    
    % For the dihedral group D_d, rho must account for 2*d elements
    if length(rho) ~= 2*d
        error('For the dihedral group D_d, the distribution rho must be of length 2*d.');
    end
    
    X = zeros(d, M);
    
    % Sample the group elements based on the distribution rho
    rho_vec = cumsum(rho);
    g_ind = rand(M, 1);
    group_elements = discretize(g_ind, [0; rho_vec]);
    
    % Calculate empirical distribution (using 1:(2*d+1) as bin edges for histcounts)
    rho_emp = histcounts(group_elements, 1:(2*d+1)); 
    rho_emp = (rho_emp / sum(rho_emp))';
    
    for m = 1 : M
        g = group_elements(m);
        
        if g <= d
            % Pure rotation (cyclic shift)
            shift_val = g - 1;
            X(:, m) = circshift(x, shift_val);
        else
            % Reflection followed by rotation
            shift_val = (g - d) - 1;
            
            % Reflect: (s*x)[l] = x[-l mod d]
            % In 1-based MATLAB indexing:
            x_reflected = [x(1); flip(x(2:end))]; 
            
            % Shift the reflected signal
            X(:, m) = circshift(x_reflected, shift_val);
        end
    end
    
    if ~exist('noisetype', 'var') || isempty(noisetype)
        noisetype = 'Gaussian';
    end
    
    switch lower(noisetype)
        case 'gaussian'
            if isreal(x)
                X = X + sigma*randn(d, M);
            else
                X = X + sigma*(randn(d, M) + 1i*randn(d, M))/sqrt(2);
            end
        case 'uniform'
            if isreal(x)
                X = X + sigma*(rand(d, M)-.5)*sqrt(12);
            else
                error('Uniform complex noise not supported yet.');
            end
        otherwise
            error('Noise type can be ''Gaussian'' or ''uniform''.');
    end
end