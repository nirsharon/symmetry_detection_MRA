function [snr1, M1_relErr, M2_relErr] = plot_pattern(n, m, N, sigma, savit)
% Zn and Zm
% N - number of observations
% sigma -- noise variance

close all;

if nargin<5
    savit = 0;
end

% perlim
circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';  % Circular mat.
F    = fft(eye(n))/sqrt(n);                     % Fourier

% the distribution
rho = rand(n,1); 
rho = rho/sum(rho);

% the signal
x     = randn(n,1);
sym_x = x;
%ratio = n/m;
for t=m:m:(n-1)
    sym_x = sym_x + circshift(x,t);
end
noise_term = sigma*randn(n, 1);
snr1 = norm(sym_x)^2/norm(noise_term)^2;
%s1 = 10^(snr(sym_x,noise_term) /10 )

% ground truth
M1 = circ(sym_x)*rho;
M2 = circ(sym_x)*diag(rho)*circ(sym_x)';

% gets the empiricals
[X, ~, ~] = generate_observations(sym_x, N, sigma, rho);

% empirical moments
emp_M1 = sum(X,2)/N;
emp_M2 = X*X'/N-sigma^2*eye(n);

M1_relErr = norm(emp_M1 - M1)/norm(M1);
M2_relErr = norm(emp_M2 - M2,'fro')/norm(M2,'fro');

% pattern ploting:

%[(F*sym_x).*(F*rho), F*M1]
%imagesc(abs([F*circ(sym_x)*diag(rho)*circ(sym_x)'*F', F*M2*F']))
figure;
stem(abs(F*emp_M1), 'filled')
if savit
    nameit = ['M1_n_',num2str(n),'_m_',num2str(m)];
    saveas(gcf,nameit,'fig');        
    saveas(gcf,nameit,'jpg');
end

figure;
imagesc(abs(F*emp_M2*F'));
colormap('pink');
colorbar

if savit
    nameit = ['M2_n_',num2str(n),'_m_',num2str(m)];
    saveas(gcf,nameit,'fig');        
    saveas(gcf,nameit,'jpg');
end


end

