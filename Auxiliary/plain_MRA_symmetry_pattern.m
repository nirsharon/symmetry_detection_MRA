% script name: "plain_MRA_symmetry_pattern"

clear;

% perlim
circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';


% parameters

n = 12;
m = 4;
N = 10000;                  % number of observations
rho   = rand(n,1); 
rho   = rho/sum(rho);
sigma = 0.1;
x     = randn(n,1);


sym_x = x;
ratio = n/m;
for t=m:m:(n-1)
    sym_x = sym_x + circshift(x,t);
end
[X, shifts,rho_emp] = generate_observations(sym_x, N, sigma, rho);
M1 = sum(X,2)/N;
M2 = X*X'/N-sigma^2*eye(n);
F = fft(eye(n))/sqrt(n);
[circ(sym_x)*rho, M1]

norm(circ(sym_x)*rho- M1)/norm(M1)

[(F*sym_x).*(F*rho), F*M1]

[circ(sym_x)*diag(rho)*circ(sym_x)', M2]

norm(circ(sym_x)*diag(rho)*circ(sym_x)'- M2,'fro')/norm(M2,'fro')

imagesc(abs([F*circ(sym_x)*diag(rho)*circ(sym_x)'*F', F*M2*F']))

