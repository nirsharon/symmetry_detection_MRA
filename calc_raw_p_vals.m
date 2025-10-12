
function[nodes_info] = calc_raw_p_vals(nodes_info, power_spectrum, sigma, N)
% This function calculates the raw p_values
% input
% power_spectrum  -- the diagonal of the second moment (power spectrum)
% sigma -- the variance
% N     -- # of observations

d            = nodes_info(1).order;
num_of_nodes = length(nodes_info);

gamma = (sigma^2*d)/N;

bdf   = N*d;                    % the basic degre of freedom
% t_Nd = chi2inv(1 - alpha, df);

% p = p = 1 - chi2cdf(chi_stat, df); % right-tailed p-value
% p = chi2cdf(chi_stat, df);         % left-tailed p-value

% traversing the graphs
for j=1:num_of_nodes

    % m is the divisor, calculate T_m
    m = nodes_info(j).order;

    % get the value of gamma for m
    %gamma_m = (1-1/m)*gamma;
    %if gamma_m==0
    %    p=0;
    %else
    Idm = calc_Idm(d, m);
    T_m = sum(power_spectrum(Idm==0));

    % should be d*N
    df = bdf*(1-1/m);
    p = 1 - chi2cdf(T_m/gamma, df); % right-tailed p-value
    
    %end
    nodes_info(j).p_value = p;
end