% script name: "robustness_test_dihedral"
% cleaning the table
close all;
clear;

% basic parameters for collapsed regime (d/m is odd)
d = 81;  
m = 9;   % d/m = 9 (odd)

tic
%--------- main parameters ----------
to_save = 1;
repeats = 20;
num_of_levels = 25;
min_sigma = 0.5;
max_sigma = 15;
%------------------------------------

% 1. Create the D_m symmetric signal
base_x = randn(d,1);
x = zeros(d, 1);
% Symmetrize over D_m (shifts by multiples of m, and their reflections)
for k = 0:(d/m - 1)
    shift_val = k * m;
    % Shift
    x_shifted = circshift(base_x, shift_val);
    % Reflect (s*x)[l] = x[-l mod d]
    x_reflected = [x_shifted(1); flip(x_shifted(2:end))]; 
    x = x + x_shifted + x_reflected;
end

% Force the signal to have exactly an average energy of 1 per entry
x = x / norm(x) * sqrt(d);

% construct the Hasse diagram graph for Dihedral
% (Assuming you have this function based on the files you listed earlier)
[G, nodes_info] = hasse_diagram_dihedral_group(d); 

% the distribution (now over 2d elements for D_d)
rho = abs(rand(2*d,1));
rho = rho/sum(rho);

M = 10000;       % # of observations
noisetype = 'gaussian';
sigma2_SNR = norm(x)^2/d;
sigma_levels  = linspace(min_sigma, max_sigma, num_of_levels);

% arrays initialization
sucsess_arr = zeros(num_of_levels,repeats);
snr_arr     = zeros(num_of_levels,1);
ps_s_succs  = zeros(num_of_levels,1);

for noise_lev = 1:num_of_levels
    sigma = sigma_levels(noise_lev);   % noise variance
    current_snr = sigma2_SNR/(sigma^2);
    snr_arr(noise_lev) = current_snr;
    
    % start testing
    for t = 1:repeats
        % make observations (needs to generate from D_d)
        % Ensure generate_observations handles 2*d states (d shifts + d shift-reflections)
        [X, shifts, rho_emp] = generate_observations_Dd(x, M, sigma, rho, noisetype);
        
        % Move to Fourier domain for moment analysis
        FX = fft(X); % d x M matrix
        
        % -------------------------------------------------------------
        % Evaluate moments for the specific hypothesis D_m
        % In practice, calc_raw_p_vals_dihedral will loop over all hypothesis nodes.
        % Here we show the math for the specific target D_m to illustrate the Wald test.
        % -------------------------------------------------------------
        
        % Step A: 2nd Order features (S_l)
        % For generic cyclic check (l not multiple of m)
        S_l_j = abs(FX).^2 - d*sigma^2; % d x M
        S_l = mean(S_l_j, 2);           % d x 1
        
        % Step B: 4th Order features (H_t)
        % J_{m,d} tuples: t = (0, m, k*m) for 1 <= k <= floor((d/m - 1)/2)
        k_max = floor(((d/m) - 1) / 2);
        q_m_dihedral = k_max;
        
        H_t = zeros(q_m_dihedral, 1);
        D_t_j_all = zeros(q_m_dihedral, M); % store for covariance calculation
        
        for k = 1:k_max
            l1 = 1;         % index 0 in MATLAB (1-based)
            l2 = m + 1;     % index m
            l3 = k*m + 1;   % index km
            
            l4_idx = mod(0 + m + k*m, d); 
            l4 = l4_idx + 1; 
            
            % d_t,j calculation
            term1 = FX(l1, :) .* FX(l2, :) .* FX(l3, :) .* conj(FX(l4, :));
            term2 = conj(FX(l1, :)) .* FX(l2, :) .* FX(l3, :) .* FX(l4, :);
            d_t_j = real(term1 + term2); % size 1 x M
            D_t_j_all(k, :) = d_t_j;
            
            D_t = mean(d_t_j);
            
            % H_t calculation
            % l4 folded representative l4_star = min(l, d-l)
            l4_star = min(l4_idx, d - l4_idx) + 1;
            
            H_t(k) = 4 * S_l(l1) * S_l(l2) * S_l(l3) * S_l(l4_star) - D_t^2;
        end
        
        % -------------------------------------------------------------
        % Step C: Assign initial p-values for all nodes
        % -------------------------------------------------------------
        % You will need to wrap the above logic into `calc_raw_p_vals_dihedral`.
        % For cyclic nodes C_m: use the Chi-squared test on S_l as you did before.
        % For dihedral nodes D_m: Construct the g_m vector and Omega_g covariance
        % matrix from S_l_j and D_t_j, invert it, and apply chi2cdf on the Wald Statistic W_m.
        
        new_nodes_info = calc_raw_p_vals_dihedral(nodes_info, FX, sigma, M);
        
        % the reject algorithm (step 5)
        [new_D] = sequential_rejection_DAG(G, new_nodes_info);
        
        % concluding (step 6)
        ind = find(new_D==0);
        if isempty(ind)
            m_hat = 1;
            symmetry_type = 'C'; 
        else
            m_hat = new_nodes_info(ind(1)).order;
            symmetry_type = new_nodes_info(ind(1)).type; % 'C' or 'D'
        end
        
        % success flag (exact match of group order and dihedral type)
        flag = 0;
        if (m_hat == m) && strcmp(symmetry_type, 'D')
            sucsess_arr(noise_lev,t) = sucsess_arr(noise_lev,t)+1;
            flag = 1;
        end
    end
    ps_s_succs(noise_lev) = flag;
end

% concluding
sucsess = sum(sucsess_arr,2)/repeats;

figure;
semilogx(snr_arr,sucsess,'LineWidth',1.5);
xlim([min(snr_arr), max(snr_arr)])
grid on;
ax = gca; % get current axes
ax.FontSize = 18;
ylabel('Success Rate')
xlabel('SNR')
title('Dihedral Symmetry Detection ($D_9$)')
set(findall(ax, '-property', 'Interpreter'), 'Interpreter', 'latex')

if to_save
    folder_name = ['success_rate_dihedral_' datestr(today,'yyyy-mm-dd')]
    if ~isfolder(folder_name)
        mkdir(folder_name)
    end
    cd(folder_name)
    name_it = ['d_', num2str(d),'_Dm_', num2str(m) ] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    exportgraphics(gcf, [name_it,'.pdf'], 'ContentType', 'vector');
    save(['data_for_',name_it]);
    cd '../'
end
toc()
% =========================================================================
% AUTOMATIC FEATURE VISUALIZATION 
% Generates a snapshot of the 2nd and 4th order features for the paper
% =========================================================================
disp('Generating feature snapshots for visual comparison...')

% Pick the highest and lowest SNR levels tested in the main loop
sigmas_to_plot = [sigma_levels(1), sigma_levels(end)];
titles = {'High SNR (Success Case)', 'Low SNR (Failure Case)'};

% Initialize the figure and the tiled layout
% 'TileSpacing', 'compact' pulls the columns and rows closer together. 
% You can also try 'none' if 'compact' still leaves too much space.
figure;
%tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
t = tiledlayout(2, 2);

% --- SPACE SETTINGS ---
t.TileSpacing = 'normal'; % Try 'normal', 'compact', or 'none'
t.Padding = 'compact';     % This controls the outer border margin
% --------------------------------

%figure('Position', [100, 100, 1200, 500]);

for i = 1:2
    sigma_val = sigmas_to_plot(i);
    snr_val = sigma2_SNR / (sigma_val^2);
    
    % Generate one fresh set of observations for this specific noise level
    [X_vis, ~, ~] = generate_observations_Dd(x, M, sigma_val, rho, noisetype);
    FX_vis = fft(X_vis);
    
    % --- Evaluate 2nd Order Features (S_l) ---
    S_l_j_vis = abs(FX_vis).^2 - d*sigma_val^2; 
    S_l_hat_vis = mean(S_l_j_vis, 2); 
    
    % --- Evaluate 4th Order Features (H_t) ---
    k_max_vis = floor(((d/m) - 1) / 2);
    H_t_hat_vis = zeros(k_max_vis, 1);
    
    for k = 1:k_max_vis
        l1 = 1; l2 = m + 1; l3 = k*m + 1;
        l4_idx = mod(0 + m + k*m, d);
        l4 = l4_idx + 1;
        
        term1 = FX_vis(l1, :) .* FX_vis(l2, :) .* FX_vis(l3, :) .* conj(FX_vis(l4, :));
        term2 = conj(FX_vis(l1, :)) .* FX_vis(l2, :) .* FX_vis(l3, :) .* FX_vis(l4, :);
        d_t_j_vis = real(term1 + term2); 
        D_t_hat_vis = mean(d_t_j_vis);
        
        l4_star_idx = min(l4_idx, d - l4_idx);
        l4_star = l4_star_idx + 1;
        
        H_t_hat_vis(k) = 4 * S_l_hat_vis(l1) * S_l_hat_vis(l2) * S_l_hat_vis(l3) * S_l_hat_vis(l4_star) - D_t_hat_vis^2;
    end
    
    % --- Plotting ---
    % Plot Power Spectrum
    nexttile(i);  % Or subplot(2, 2, i) if you reverted to the old method
    stem(0:(d-1), S_l_hat_vis, 'filled', 'MarkerSize', 4);
    
    % FIX 1: Added 'FontWeight', 'normal' to remove the bolding
    title(sprintf('%s\nPower Spectrum (SNR = %.3f)', titles{i}, snr_val), 'FontWeight', 'normal');
    
    xlabel('Frequency Index (l)');
    ylabel('$\hat{S}_l$', 'Interpreter', 'latex', 'FontSize', 14); 
    xlim([0 d-1]);
    grid on;
    
    % Plot 4th Order Features
    nexttile(i + 2); 
    stem(1:k_max_vis, H_t_hat_vis, 'filled', 'MarkerSize', 6, 'Color', 'r');
    
    % FIX: The empty string goes SECOND to push the title UP
    title({'4th-Order Dihedral Features ($\hat{H}_t$)', ' '}, 'Interpreter', 'latex', 'FontSize', 14);
    
    xlabel('Tuple Index (k)');
    ylabel('$\hat{H}_t$', 'Interpreter', 'latex', 'FontSize', 14); 
    xlim([0 k_max_vis+1]);
    grid on;
end

% Optional: Save this new figure too
if to_save
    cd(folder_name)
    saveas(gcf, [name_it, '_features'] ,'fig');
    %saveas(gcf,[name_it, '_features'] ,'jpg');
    exportgraphics(gcf, [name_it,'_features.pdf'], 'ContentType', 'vector');
    cd '../'
end