function [nodes_info] = calc_raw_p_vals_dihedral(nodes_info, FX, sigma, N)
% This function calculates the raw p_values for both Cyclic and Dihedral DAG nodes
% inputs:
% nodes_info -- struct containing the DAG nodes (must have 'order' and 'type' fields)
% FX         -- d x N matrix of Fourier transformed observations (fft(X))
% sigma      -- the noise standard deviation
% N          -- # of observations

d = size(FX, 1);
num_of_nodes = length(nodes_info);
gamma = (sigma^2 * d) / N;
bdf = N * d; % basic degree of freedom for cyclic chi-squared test

% Precompute the power spectrum for the Cyclic tests
power_spectrum = mean(abs(FX).^2, 2);

% Precompute the second-order sample features (S_l,j) for Dihedral tests
S_l_j = abs(FX).^2 - d * sigma^2; % d x N
S_l_hat = mean(S_l_j, 2);         % d x 1

% Traverse the graph
for j = 1:num_of_nodes
    m = nodes_info(j).order;
    type = nodes_info(j).type; % Assumes nodes have a type field: 'C' or 'D'

    if strcmp(type, 'C')
        % ---------------------------------------------------------
        % CYCLIC NODE TESTING (Original Logic)
        % ---------------------------------------------------------
        Idm = calc_Idm(d, m); % Your helper function mapping multiples of m
        T_m = sum(power_spectrum(Idm == 0));

        df = bdf * (1 - 1/m);
        p = 1 - chi2cdf(T_m / gamma, df); % right-tailed p-value
        nodes_info(j).p_value = p;

    elseif strcmp(type, 'D')
        % ---------------------------------------------------------
        % DIHEDRAL NODE TESTING (Collapsed Regime, Section 5.3.4)
        % ---------------------------------------------------------

        % 1. Define index sets (0-based math converted to 1-based MATLAB)
        % I_{m,d}^{c,+} : Indices not multiple of m, from 0 to floor(d/2)
        I_md_c_plus = [];
        for l_idx = 0:floor(d/2)
            if mod(l_idx, m) ~= 0
                I_md_c_plus = [I_md_c_plus; l_idx + 1];
            end
        end
        % Absolute guarantee it is a column vector
        I_md_c_plus = I_md_c_plus(:);

        % 2. Evaluate Cyclic Part of g_m
        g_m_cyclic = S_l_hat(I_md_c_plus);
        g_m_cyclic = g_m_cyclic(:); % Absolute guarantee it is a column vector

        % 3. Evaluate Dihedral Part of g_m (J_{m,d} tuples)
        k_max = floor(((d/m) - 1) / 2);
        g_m_dihedral = zeros(k_max, 1);

        % Pre-allocate D_t_j for covariance building
        D_t_j_all = zeros(k_max, N);
        D_t_hat_all = zeros(k_max, 1);

        for k = 1:k_max
            l1 = 1;               % 0 mod d
            l2 = m + 1;           % m mod d
            l3 = k*m + 1;         % km mod d

            l4_idx = mod(0 + m + k*m, d);
            l4 = l4_idx + 1;

            % Compute d_{t,j} across all N observations
            term1 = FX(l1, :) .* FX(l2, :) .* FX(l3, :) .* conj(FX(l4, :));
            % Correctly evaluate 2 * Real(term1) to form the invariant D_t
            d_t_j = 2 * real(term1); 
            D_t_j_all(k, :) = d_t_j;

            D_t_hat = mean(d_t_j);
            D_t_hat_all(k) = D_t_hat;

            % Folded representative for S_l (Eq 32)
            l4_star_idx = min(l4_idx, d - l4_idx);
            l4_star = l4_star_idx + 1;

            % Compute H_t (Eq 34 sample version)
            g_m_dihedral(k) = 4 * S_l_hat(l1) * S_l_hat(l2) * S_l_hat(l3) * S_l_hat(l4_star) - D_t_hat^2;
        end

        % Combine into full test vector
        g_m = [g_m_cyclic; g_m_dihedral];
        q_m = length(g_m); % Degrees of freedom for Wald test

        % 4. Build empirical covariance matrix \hat{\Omega}_g (Eq 41)
        % Build the \hat{\psi}_j matrix (q_m x N)
        psi_hat = zeros(q_m, N);

        % Cyclic rows - safely handled with bsxfun
        num_cyclic = length(I_md_c_plus);
        if num_cyclic > 0
            S_l_j_sub = S_l_j(I_md_c_plus, :);
            % bsxfun automatically handles the dimension matching for subtraction
            psi_hat(1:num_cyclic, :) = bsxfun(@minus, S_l_j_sub, g_m_cyclic);
        end

        % Dihedral rows (Eq 40 implementation)
        for k = 1:k_max
            l1 = 1;
            l2 = m + 1;
            l3 = k*m + 1;
            l4_idx = mod(0 + m + k*m, d);
            l4_star = min(l4_idx, d - l4_idx) + 1;

            % Delta method expansion parts
            part_0 = S_l_hat(l2) * S_l_hat(l3) * S_l_hat(l4_star) * (S_l_j(l1, :) - S_l_hat(l1));
            part_m = S_l_hat(l1) * S_l_hat(l3) * S_l_hat(l4_star) * (S_l_j(l2, :) - S_l_hat(l2));
            part_km = S_l_hat(l1) * S_l_hat(l2) * S_l_hat(l4_star) * (S_l_j(l3, :) - S_l_hat(l3));
            part_4star = S_l_hat(l1) * S_l_hat(l2) * S_l_hat(l3) * (S_l_j(l4_star, :) - S_l_hat(l4_star));

            S_derivs = 4 * (part_0 + part_m + part_km + part_4star);
            D_derivs = -2 * D_t_hat_all(k) * (D_t_j_all(k, :) - D_t_hat_all(k));

            psi_hat(num_cyclic + k, :) = S_derivs + D_derivs;
        end

        % --- NEW NORMALIZATION FIX ---
        % Calculate the standard deviation of each feature row (add tiny epsilon to prevent /0)
        scale_factors = sqrt(sum(psi_hat.^2, 2) / N) + 1e-10; 
        
        % Scale g_m and psi_hat so every feature has a variance of ~1
        g_m_scaled = g_m ./ scale_factors;
        psi_hat_scaled = bsxfun(@rdivide, psi_hat, scale_factors);
        
        % Calculate standardized sample covariance matrix
        Omega_g_hat_scaled = (1/N) * (psi_hat_scaled * psi_hat_scaled');
        
        % --- NEW RIDGE REGULARIZATION --- OLD!
        % Add a tiny diagonal ridge to prevent highly collinear 4th-order features 
        % from creating near-zero eigenvalues that explode during inversion.
        % Omega_g_hat_scaled = Omega_g_hat_scaled + 1e-4 * eye(q_m);
        
        % --- DELICATE SVD TRUNCATION ---
        % 1. Define a singular value tolerance threshold. 
        % You can tune this. 

        % Set tolerance relative to the largest singular value of the matrix
        % This drops any component that is 1 million times smaller than the primary variance
        % --- DELICATE SVD TRUNCATION ---
        svd_tol = 1e-6 * norm(Omega_g_hat_scaled); 
        
        % 1. Dynamically compute the effective degrees of freedom
        effective_dof = rank(Omega_g_hat_scaled, svd_tol);
        
        % 2. Compute the stable Wald Statistic
        W_m = N * g_m_scaled' * pinv(Omega_g_hat_scaled, svd_tol) * g_m_scaled;
        
        % 6. Calculate p-value dynamically using effective rank
        if effective_dof > 0
            p = chi2cdf(W_m, effective_dof, 'upper'); 
        else
            p = 0; 
        end
        nodes_info(j).p_value = p;
    end

end
end