% script name: "robustness_test"


% cleaning the table
close all;
clear;

% basic parameters
d = 84;  % 840
m = 21;  % 56

tic

%--------- main parameters ----------
to_save = 1;
repeats = 20;
num_of_levels = 25;

min_sigma = 1;
max_sigma = 50;
%------------------------------------

% the signal
base_x = randn(d,1);
x = get_symmetric_Zm_signal_Zn(d, d/m, base_x);

% construct the Hasse diagran graph
[G, nodes_info] = hasse_diagram_cyclic_group(d);

% the distribution
rho = abs(rand(d,1));
rho = rho/sum(rho);

M = 10000;       % # of observations
noisetype = 'gaussian';

sigma2_SNR = norm(x)^2/d;
sigma_levels  = linspace(min_sigma, max_sigma, num_of_levels);

% arrays initialization
sucsess_arr = zeros(num_of_levels,repeats);
snr_arr     = zeros(num_of_levels,1);
ps_samples  = zeros(d,num_of_levels);
ps_s_succs  = zeros(num_of_levels,1);

for noise_lev = 1:num_of_levels
    sigma = sigma_levels(noise_lev);   % noise variance
    current_snr = sigma2_SNR/(sigma^2);
    snr_arr(noise_lev) = current_snr;

    % start testing
    for t = 1:repeats

        % make observations
        [X, shifts, rho_emp] = generate_observations_Zn(x, M, sigma, rho, noisetype);

        % the empirical moments
        mu1 = mean(X,2);
        mu2 = (1/M)*(X*X');%-sigma^2*eye(n);
        F = fft(eye(d));
        ps = real(diag(F*mu2*F'));

        % assign initial p_values
        new_nodes_info = calc_raw_p_vals(nodes_info, ps, sigma, M);

        % the reject algorithm (step 5)
        [new_D] = sequential_rejection_DAG(G, new_nodes_info);

        % concluding (step 6)
        ind = find(new_D==0);
        if isempty(ind)
            m_hat = 1;
           % fprintf('There is no symmetry presented in this signal .\n')
        else
            m_hat = new_nodes_info(ind(1)).order;
          %  fprintf('The symmetry presented in this signal is m=%d.\n', m_hat);
        end

        % sucess
        flag = 0;
        if m_hat==m
            sucsess_arr(noise_lev,t) = sucsess_arr(noise_lev,t)+1;
            flag = 1;
        end
    end
    ps_samples(:,noise_lev) = ps;
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
set(findall(ax, '-property', 'Interpreter'), 'Interpreter', 'latex')

if to_save
    folder_name = 'success_rate';
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end

    name_it = ['d_', num2str(d),'_m_', num2str(m) ] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    exportgraphics(gcf, [name_it,'.pdf'], 'ContentType', 'vector');
    save(['data_for_',name_it]);
    cd '../'
end

plot_power_spectrum

toc()


