%% run_fwer_distfree.m
% Two aspects: (1) empirical DAG-level FWER control, and
%              (2) distribution-free behaviour: uniform vs strongly non-uniform rho.
%
% FAMILY-WISE ERROR = the procedure rejects AT LEAST ONE hypothesis that is
% actually a true invariance of the ground-truth signal. Valid strong control
% keeps this <= alpha for every SNR and every rho (typically strictly below).

close all; clear; rng(1);
tic
%% ---------------- parameters ----------------
to_save   = 1;
d         = 81;                          % odd => collapsed regime
alpha0    = 0.05;
T         = 500;                         % Monte-Carlo trials per point
N         = 2e4;                         % observations per trial
noisetype = 'gaussian';
num_snr   = 10;
snr_grid  = logspace(-1, 2, num_snr);    % SNR := (||x||^2/d)/sigma^2, with ||x||^2 = d

configs   = { struct('type','C','m',9), struct('type','D','m',9) };  % Sym(x) = C_9, D_9
cfg_name  = {'true C_9','true D_9'};
rho_types = {'uniform','nonuniform'};

[G, nodes_template] = hasse_diagram_dihedral_group(d);
styles = {'-o','--s'};
se = @(p) sqrt(max(p.*(1-p),0)/T);       % binomial Monte-Carlo SE

only2 = 1;
if ~only2
%% ============ Experiment 1: empirical FWER vs SNR ============
FWER = nan(num_snr, numel(rho_types), numel(configs));
POW  = nan(num_snr, numel(rho_types), numel(configs));   % correct-Sym id (side info)

for c = 1:numel(configs)
    cfg  = configs{c};
    x    = make_symmetric_signal(d, cfg.m, cfg.type);
    tn   = true_null_mask(nodes_template, cfg.type, cfg.m);
    sig2 = norm(x)^2 / d;
    for ir = 1:numel(rho_types)
        rho = make_rho(rho_types{ir}, d);
        for is = 1:num_snr
            sigma = sqrt(sig2 / snr_grid(is));
            fw = 0; pw = 0;
            parfor t = 1:T
                nodes = compute_pvalues(x, N, sigma, rho, noisetype, nodes_template);
                [rej, mhat, that] = decide(G, nodes, alpha0);
                if any(rej(tn)==1), fw = fw + 1; end
                if strcmp(that,cfg.type) && mhat==cfg.m, pw = pw + 1; end
            end
            FWER(is,ir,c) = fw/T;  POW(is,ir,c) = pw/T;
            fprintf('cfg %-8s | rho %-10s | SNR %7.3f | FWER %.3f | id-rate %.3f\n', ...
                cfg_name{c}, rho_types{ir}, snr_grid(is), FWER(is,ir,c), POW(is,ir,c));
        end
    end
end

figure('Name','FWER vs SNR'); hold on;
for c = 1:numel(configs)
    for ir = 1:numel(rho_types)
        p = FWER(:,ir,c);
        errorbar(snr_grid, p, se(p), styles{ir}, 'LineWidth',1.3, ...
            'DisplayName', sprintf('%s, %s \\rho', cfg_name{c}, rho_types{ir}));
    end
end
plot(snr_grid([1 end]), [alpha0 alpha0], 'k:', 'LineWidth',1.2, 'HandleVisibility','off');
text(snr_grid(1), alpha0, '  \alpha = 0.05', 'VerticalAlignment','bottom');
set(gca,'XScale','log','FontSize',14); grid on;
xlabel('SNR'); ylabel('empirical FWER'); ylim([0, max(0.12,1.4*alpha0)]);
legend('Location','best','Interpreter','tex');
%title('DAG-level FWER: empirical FWER across SNR and \rho','Interpreter','tex');
end
%% ============ Experiment 2: FWER calibration vs nominal alpha ============
snr_cal    = 10;
alpha_grid = linspace(0.01, 0.20, 12);
cfg = configs{1};                                   % cyclic config
x   = make_symmetric_signal(d, cfg.m, cfg.type);
tn  = true_null_mask(nodes_template, cfg.type, cfg.m);
sigma = sqrt((norm(x)^2/d) / snr_cal);

FWER_cal = nan(numel(alpha_grid), numel(rho_types));
for ir = 1:numel(rho_types)
    rho = make_rho(rho_types{ir}, d);
    err = zeros(numel(alpha_grid),1);
    parfor t = 1:T
        RandStream.setGlobalStream(RandStream('Threefry','Seed',1000*ir+t));
        nodes = compute_pvalues(x, N, sigma, rho, noisetype, nodes_template);  % p-values once
        e = zeros(numel(alpha_grid),1);
         for ia = 1:numel(alpha_grid)
             rej = decide(G, nodes, alpha_grid(ia));
             e(ia) = any(rej(tn)==1);
         end
    err = err + e;        % vector + reduction, valid
    end
    FWER_cal(:,ir) = err/T;
end

figure('Name','FWER calibration'); hold on;
plot(alpha_grid, alpha_grid, 'k:', 'LineWidth',1.2, 'DisplayName','nominal (y = x)');
for ir = 1:numel(rho_types)
    p = FWER_cal(:,ir);
    errorbar(alpha_grid, p, se(p), styles{ir}, 'LineWidth',1.3, ...
        'DisplayName', sprintf('empirical, %s \\rho', rho_types{ir}));
end
set(gca,'FontSize',14); grid on; axis square;
xlabel('nominal \alpha'); ylabel('empirical FWER');
title(sprintf('FWER calibration (%s, SNR = %g)', cfg_name{1}, snr_cal),'Interpreter','tex');
legend('Location','best','Interpreter','tex');

%% ---------------- save ----------------
if to_save
    folder = ['fwer_distfree_v3' datestr(now,'yyyy-mm-dd')];
    if ~isfolder(folder), mkdir(folder); end
    exportgraphics(figure(1), 'fwer_vs_snr.jpg');
    exportgraphics(figure(2), 'fwer_calibration.jpg');
    saveas(figure(1), fullfile(folder,'fwer_vs_snr'), 'fig');
    exportgraphics(figure(1), fullfile(folder,'fwer_vs_snr.pdf'), 'ContentType','vector');
    saveas(figure(2), fullfile(folder,'fwer_calibration'), 'fig');
    exportgraphics(figure(2), fullfile(folder,'fwer_calibration.pdf'), 'ContentType','vector');
    save(fullfile(folder,'fwer_distfree_data.mat'), ...
        'FWER','POW','FWER_cal','snr_grid','alpha_grid','d','T','N','configs','rho_types');
end
toc()
%% ==================== local functions ====================
function x = make_symmetric_signal(d, m, type)
% Exact symmetry C_m ('C') or D_m ('D'). Convention: C_m = <r^{d/m}>, so
% C_m-invariant <=> period d/m <=> Fourier support on multiples of m; hence
% average over shifts that are multiples of (d/m). NB: v2 shifts by multiples
% of m, which agrees only when m = d/m (e.g. d=81,m=9). Use d/m in general.
    assert(mod(d,m)==0,'m must divide d');
    step = d/m; base = randn(d,1); x = zeros(d,1);
    for k = 0:m-1, x = x + circshift(base, k*step); end
    if strcmp(type,'D'), x = x + [x(1); flipud(x(2:end))]; end   % reflection s (j=0)
    x = x / norm(x) * sqrt(d);
end

function rho = make_rho(type, d)
% Probability vector over the 2d dihedral elements (order matching generate_observations_Dd).
    switch type
        case 'uniform',    rho = ones(2*d,1)/(2*d);
        case 'nonuniform', v = rand(2*d,1).^8; v(randi(2*d)) = v(randi(2*d)) + 5*max(v); rho = v/sum(v);
        otherwise, error('unknown rho type');
    end
end

function tn = true_null_mask(nodes, star_type, star_m)
% H=(type,m) is a true null iff H<=star:  C_m<=C_{m*} or D_{m*} iff m|m* ;
% D_m<=D_{m*} iff m|m* ; D_m<=C_{m*} never.
    n = numel(nodes); tn = false(n,1);
    for j = 1:n
        divides = (mod(star_m, nodes(j).order)==0);
        switch nodes(j).type
            case 'C', tn(j) = divides;
            case 'D', tn(j) = divides && strcmp(star_type,'D');
        end
    end
end

function nodes = compute_pvalues(x, N, sigma, rho, noisetype, nodes_template)
    [X,~,~] = generate_observations_Dd(x, N, sigma, rho, noisetype);
    nodes = calc_raw_p_vals_dihedral(nodes_template, fft(X), sigma, N);
end

function [rej, mhat, that] = decide(G, nodes, alpha)
% Sets per-node level and runs the DAG procedure. ASSUMES sequential_rejection_DAG
% reads each node's .alpha field (as hasse_diagram_dihedral_group sets it). If it
% instead takes a global alpha, change to sequential_rejection_DAG(G, nodes, alpha).
    [nodes.alpha] = deal(alpha);
    rej = logical(sequential_rejection_DAG(G, nodes)); rej = rej(:);
    ind = find(rej==0);
    if isempty(ind), mhat = 1; that = 'C';
    else, mhat = nodes(ind(1)).order; that = nodes(ind(1)).type; end
end