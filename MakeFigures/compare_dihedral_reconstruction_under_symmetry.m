% script name: "compare_dihedral_reconstruction_under_symmetry"
%
% Apply MRA-DIHEDRAL-EM for special signals with symmetry
%
% NS, July 2025

clear
clc;
close all;
tic

%% main parameters
d = 60;    % signal length
m = 20; j = 2;   % signal's dihedral symmetry
N = 10^4;        % number of observations

repeat_trial = 5;   % number of repeated tests
to_save = 1;        % saving the outputs of the script
to_plot_signal = 1;

%% EM parameters
tol      = 1e-3;          % absolute fitting tolerance
max_iter = 200;           % maximum number of iterations

%% Preparing the Signals
t = linspace(0, 2*pi, d);
%signal = sin(t) + 0.5*sin(3*t + pi/4) + 0.25*sin(5*t + pi/2) + 0.01*randn(1, d);
signal = cos(1.5*t + pi/3) + 0.6*sin(4*t + pi/5) + 0.3*cos(2*t + pi/7) + 0.01*randn(1,d);
signal = signal(:);

sigma2_SNR   = norm(signal)^2/d;
signal_c_sym = get_symmetric_Zm_signal_Dn(d, m, signal);
signal_D_sym = get_symmetric_Dmj_signal_Dn(d, m, j, signal);

% Plotting
if to_plot_signal
    dom = 0:(d-1);
    figure;
    hold on
    plot(dom, signal, '-o', 'LineWidth', 1.5);
    plot(dom, signal_c_sym, 'LineWidth', 1.8);
    plot(dom, signal_D_sym, 'LineWidth', 1.8);
    name1 = ['{C}_{', num2str(m), '} symmetry'];
    %name2 = ['{D}_{', num2str(m),',',num2str(j), '} symmetry'];
    name2 = ['{D}_{', num2str(m), '} symmetry'];
    h = legend('Asymmetric signal', name1 , name2,'Location','northeast');
    pos = get(h, 'Position');
    pos(1) = pos(1) + 0.05; % adjust this value to your preference
    set(h, 'Position', pos);

    %xlabel('x');
    %ylabel('y');
    xlim([min(dom), max(dom)])
    ax = gca; % get current axes
    ax.FontSize = 18;
    set(findall(ax, '-property', 'Interpreter'), 'Interpreter', 'latex')
    grid on;

    if to_save
        folder_name = 'sym_reconst';
        if isfolder(folder_name)
            cd(folder_name);
        else
            mkdir(folder_name)
            cd(folder_name)
        end

        name_it = ['signals_d_', num2str(d),'_m_', num2str(m),'_j_', num2str(j) ] ;
        saveas(gcf, name_it ,'fig');
        saveas(gcf,name_it,'jpg');
        exportgraphics(gcf, [name_it,'.pdf'], 'ContentType', 'vector');
        cd '../'
    end
end

%% parallel option for parallel run
parallel_pool = 0;
if parallel_pool
    parallel_nodes = 4;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

%% setting the test

% noise and trials parameters
length_sigma = 5;           % number of dfferen sigma values
sigma_vec    = logspace(-2,.45,length_sigma);
snr_arr      = sigma2_SNR./(sigma_vec.^2);

% initial error arrays
run_em  = zeros(repeat_trial, length_sigma);
err_em_c = zeros(repeat_trial, length_sigma);
err_em_d = zeros(repeat_trial, length_sigma);

% repeat trials loop
for trial_num = 1:repeat_trial

    p = randn(2*d, 1).^2; p = p/sum(p);   % random distribution

    [err_em(trial_num,:)]   = run_dihedral_EM(signal, p, N, length_sigma, sigma_vec, max_iter, tol);
    [err_em_c(trial_num,:)] = run_dihedral_EM(signal_c_sym, p, N, length_sigma, sigma_vec, max_iter, tol);
    [err_em_d(trial_num,:)] = run_dihedral_EM(signal_D_sym, p, N, length_sigma, sigma_vec, max_iter, tol);
end

err_em = mean(err_em,1);
err_em_c = mean(err_em_c,1);
err_em_d = mean(err_em_d,1);
%% plotting

figure;
hold on;
plot(snr_arr, err_em,'LineWidth',1.6);
plot(snr_arr, err_em_c,'LineWidth',1.8);
plot(snr_arr, err_em_d,'LineWidth',1.8);

ax = gca; % get current axes
ax.FontSize = 18;
set(findall(ax, '-property', 'Interpreter'), 'Interpreter', 'latex')

name1 = ['{C}_{', num2str(m), '} symmetry'];
%name2 = ['{D}_{', num2str(m),',',num2str(j), '} symmetry'];
name2 = ['{D}_{', num2str(m), '} symmetry'];
legend('No Symmetry', name1 , name2,'Location','best');
set(ax, 'YScale', 'log')
set(ax, 'XScale', 'log')
set(ax,"XLim",[min(snr_arr),max(snr_arr)])
xlabel('SNR');
ylabel('RMSE');
grid on;

if to_save
    folder_name = 'sym_reconst';
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end

    name_it = ['compare_signals_d_', num2str(d),'_m_', num2str(m),'_j_', num2str(j) ] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    exportgraphics(gcf, [name_it,'.pdf'], 'ContentType', 'vector');
    save(['data_for_',name_it]);
    cd '../'
end

toc()

