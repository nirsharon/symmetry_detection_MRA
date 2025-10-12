% script name: plot_rejection_process
%

% cleaning the table
close all;
clear;

% basic parameters
n = 90;  % 840
m = 10;  % 56

to_save = 1;


% the signal
base_x = randn(n,1);
x = get_symmetric_Zm_signal_Zn(n, n/m, base_x);

% the distribution 
rho = abs(rand(n,1)); 
rho = rho/sum(rho);

% make observations 
M = 10000;       % # of observations
sigma = 0.5;   % noise variance
noisetype = 'gaussian';
[X, shifts, rho_emp] = generate_observations_Zn(x, M, sigma, rho, noisetype);

% the empirical moments 
mu1 = mean(X,2);
mu2 = (1/M)*(X*X')-sigma^2*eye(n);
F = fft(eye(n));
ps = real(diag(F*mu2*F'));

% construct the Hasse diagran graph 
[G, nodes_info] = hasse_diagram_cyclic_group(n);

% assign initial p_values 
alpha = 0.05;
new_nodes_info = calc_raw_p_vals(nodes_info, ps, sigma, M);

% the reject algorithm (step 5)
[new_D] = visual_sequential_rejection_DAG(G, new_nodes_info, to_save);

% concluding (step 6)
ind = find(new_D==0);
if isempty(ind)
    fprintf('There is no symmetry presented in this signal .\n')
else
    m_hat = new_nodes_info(ind(1)).order;
    fprintf('The symmetry presented in this signal is m=%d.\n', m_hat);
end
