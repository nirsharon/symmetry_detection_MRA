% script name: "run_spectral_alg"
% ==== taken from "main" =====
%%%%%%%%%%% EVEN CASE %%%%%%%%

clear
clc; 
close all; 

%% parameters
L = 10;   % signal length
N = 10^6; % number of observations
sigma = 0.01; %.01;

%% data
true_x   = randn(L, 1);    % underlying signal

true_rho = rand(2*L, 1);   % rand distribution
true_rho = true_rho/sum(true_rho); % raw stochastic
p = true_rho(1:L); q = true_rho(L+1:end);  % notational 

% imposing the assumption on the distribution
if mod(L,2)==0    % even length
    split = rand;
    q(1:2:end) = rand;
    q(2:2:end) = rand;
    q = q/(sum(q)/(1-sum(p)));
    q0 = q(1); 
    q1 = q(2); 
    
    p = (p + reverse(p))/2;
end
true_rho(L+1:end) = q;
true_rho(1:L) = p;

[X, shifts] = generate_observations(true_x, N, sigma, true_rho); 

%% sanity check

circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';
Cx = circ(true_x);

F  = 1/(sqrt(L))*fft(eye(L));
D_Fx = diag(F*true_x);

M2 = Cx*diag(p)*Cx.' + Cx.'*diag(q)*Cx;
M2_fourier = F*M2*F';

% %----------
% Fp = F*p;
% Fq = F*q;
% 
% Fp_prime = Fp;
% Fp_prime(L/2+1)=0;  
% delta_p = Fp - Fp_prime;
% 
% Fq_prime = Fq;
% Fq_prime(L/2+1)=0; 
% delta_q = Fq - Fq_prime;
% 
% % 1
% norm(M2_fourier/sqrt(L) - D_Fx*(circ(Fp))*(D_Fx') - D_Fx'*(circ(Fq))*(D_Fx))
% %2
% norm(M2_fourier/sqrt(L) - D_Fx*(circ(Fp_prime))*(D_Fx') ...
%    - D_Fx*(circ(delta_p))*(D_Fx') - D_Fx'*(circ(Fq_prime))*(D_Fx) ...
%    -D_Fx'*(circ(delta_q))*(D_Fx))
% %3
% norm(M2_fourier/sqrt(L) - D_Fx*(circ(Fp_prime))*(D_Fx') ...
%    - D_Fx*(circ(delta_p))*(D_Fx') - D_Fx'*(diag(sqrt(L)*((q0+q1)/2)))*(D_Fx) ...
%    -D_Fx'*(circ(delta_q))*(D_Fx))
% %4
% vv = zeros(L,1); vv(1+L/2)=1;
% index_mat = circ(vv);
% norm(M2_fourier/sqrt(L) - D_Fx*(circ(Fp_prime))*(D_Fx') ...
%    - D_Fx*(circ(delta_p))*(D_Fx') - D_Fx'*(diag(sqrt(L)*((q0+q1)/2)))*(D_Fx) ...
%    - D_Fx'*(index_mat*(sqrt(L)*((q0-q1)/2)))*(D_Fx) )
% %5
% norm(M2_fourier/sqrt(L) ... 
%    - (D_Fx*(circ(Fp_prime))*(D_Fx') + D_Fx'*(diag(sqrt(L)*((q0+q1)/2)))*(D_Fx)) ...
%    - (D_Fx'*(index_mat*(sqrt(L)*((q0-q1)/2)))*(D_Fx) + D_Fx*(index_mat*Fp(L/2+1))*(D_Fx') ) ...
% )
% %6
% ind_mat2 = ones(L)-index_mat;
% M2_tilde = zeros(L);
% M2_tilde(ind_mat2==1)= M2_fourier(ind_mat2==1);
% norm(M2_tilde/sqrt(L) - (D_Fx*(circ(Fp_prime))*(D_Fx') + D_Fx*(diag(sqrt(L)*((q0+q1)/2)))*(D_Fx)'))
% %7
% norm(M2_tilde/sqrt(L) - D_Fx*( circ(Fp_prime)+ diag(sqrt(L)*((q0+q1)/2))*eye(L)  )*(D_Fx)' )
% %8
% norm(F'*M2_tilde*F/sqrt(L) - F'*D_Fx*( circ(Fp_prime)+ diag(sqrt(L)*((q0+q1)/2))*eye(L)  )*(D_Fx)'*F )
% %9
% Fpp = Fp_prime;
% Fpp(1) = Fpp(1) + sqrt(L)*(q0+q1)/2;
% norm(F*M2_tilde*F'/sqrt(L) - F*D_Fx*circ(Fpp)*(D_Fx)'*F' )
% norm(F*circ(Fpp)*F'-diag(F'*Fpp)*sqrt(L))
% norm(F*(D_Fx)'*F'-Cx/sqrt(L))
% norm(sqrt(L)*F*(D_Fx)*F'-Cx')
% norm(F*M2_tilde*F' - Cx'*diag(F'*Fpp)*Cx)
% %10
% norm_mat = diag(abs(F*true_x).^(-1));
% conj_mat = F*norm_mat;
% 
% V = F'*D_Fx*norm_mat*F;
% norm(V*V'-eye(L))  % orthogonal
% M_ext = conj_mat'*M2_tilde*conj_mat;
% 
% norm(M2_tilde/sqrt(L) - D_Fx*( circ(Fp_prime)+ diag(sqrt(L)*((q0+q1)/2))*eye(L)  )*(D_Fx)' );
% V'*M_ext*V
% 
% % 11
% Fpp = Fp_prime; Fpp(1) = Fpp(1) + sqrt(L)*(q0+q1)/2;
% circ(Fpp) == (circ(Fp_prime)+ eye(L)*(sqrt(L)*(q0+q1)/2));
% norm(F*circ(Fpp)*F' - sqrt(L)*diag(F'*Fpp))
% norm(F*(D_Fx)'*F'-Cx/sqrt(L))
% 
% norm(F*(D_Fx)'*F'*F*circ(Fpp)*F'*F*(D_Fx)*F' - Cx/sqrt(L)*diag(F'*Fpp)*Cx')
% norm(F*M2_tilde*F' - Cx'*diag(F'*Fpp)*Cx)
% norm(sqrt(L)*F*(D_Fx)'*F'-Cx)
% Cx_ortho = (1/sqrt(L)*F*norm_mat*F')*Cx;
% norm(Cx_ortho*Cx_ortho'-eye(L))
% conj_mat = (1/sqrt(L)*F*norm_mat*F');
% norm(conj_mat'*Cx'*diag(F'*Fpp)*Cx*conj_mat - Cx_ortho'*diag(F'*Fpp)*Cx_ortho)
% norm(conj_mat'*F*M2_tilde*F'*conj_mat - Cx_ortho'*diag(F'*Fpp)*Cx_ortho)
% conj_M2 = F'*conj_mat;
% norm(conj_M2-1/sqrt(L)*norm_mat*F')
% 
% norm(conj_M2'*M2_tilde*conj_M2 - Cx_ortho'*diag(F'*Fpp)*Cx_ortho)
% 
% M_ext = conj_M2'*M2_tilde*conj_M2;
% V     = Cx_ortho'; % columns are eigenvectors
% %V'*M_ext*V

% % ------------------- dist calc. ignore --------
% x = F*F'*Fpp;
% Fpp = F*p; Fpp(1) = Fpp(1) + sqrt(L)*(q0+q1)/2;
% Fpp(1) == p(1) + sqrt(L)*(q0+q1)/2; 
% p(L/2) + sqrt(L)*(q0-q1)/2
% y = p+q;
% y(1) == p(1) + q1
% y(L/2) == p(L/2) + q0


%% run spectral algorithm
[est_x, ~] = spectral_alg(X, sigma, true_x, true_rho);

%% plotting 
x_est_al = align_to_reference(est_x, true_x);
err = norm(x_est_al - true_x)/norm(true_x);
fprintf('Relative error = %.4g\n', err);

figure; 
%subplot(121);
hold on; 
stem(true_x); 
stem(x_est_al);
legend({'$x$', '$\hat{x}$ (est)'},'Interpreter','latex');

% subplot(122);
% stem(true_rho); 
% stem(est_dist);
% legend('\rho', '\hat{\rho}'); 
