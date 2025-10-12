% script name: run_Zm_example
%
%==================================
% This script runs a full Zm example where we:
%
% 1. generate observations
% 2. construct moments
% 3. establish Hasse diagram (graph tree)
% 4. assign initial p-values
% 5. run iterative reject
% 6. print the conclusion
%
%==================================

% cleaning the table
close all;
clear;

% basic parameters
n = 60;
m = 10;

% the signal
base_x = randn(n,1);
x = get_symmetric_Zm_signal_Zn(n, n/m, base_x);

% the distribution 
rho = abs(rand(n,1)); 
rho = rho/sum(rho);

% make observations (step 1)
M = 10000;       % # of observations
sigma = 0.5;   % noise variance
noisetype = 'gaussian';
[X, shifts, rho_emp] = generate_observations_Zn(x, M, sigma, rho, noisetype);

% the analytical moments (just for comparison)
Cx = circulant(x);
M1 = Cx*rho;
M2 = Cx*diag(rho)*Cx';

% the empirical moments (step 2)
mu1 = mean(X,2);
mu2 = (1/M)*(X*X'); %-sigma^2*eye(n);
F = fft(eye(n));
ps = real(diag(F*mu2*F'));

% construct the Hasse diagran graph (step 3)
[G, nodes_info] = hasse_diagram_cyclic_group(n);

% assign initial p_values (step 4)
alpha = 0.05;
new_nodes_info = calc_raw_p_vals(nodes_info, ps, sigma, M);

% plot current graph 
plot_the_hasse_diagram(G, new_nodes_info)

% the reject algorithm (step 5)
 [new_D] = sequential_rejection_DAG(G, new_nodes_info);

% concluding (step 6)
ind = find(new_D==0);
if isempty(ind)
    fprintf('There is no symmetry presented in this signal .\n')
else
    m_hat = new_nodes_info(ind(1)).order;
    fprintf('The symmetry presented in this signal is m=%d.\n', m_hat);
end
