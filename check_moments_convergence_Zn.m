% script name: check_moments_convergence_Zn
%
% checking the code-generated moments against the analytic formulas

clear
clc; 
close all; 
tic

%% parameters
L = 5;   % signal length
num_of_trials = 10;
N_arr = floor(logspace(2,4.5,num_of_trials)); % observations #, ranges from 100 to 10^5
sigma = 0.001; %.01;

%% data
x  = randn(L, 1);   % underlying signal

p = rand(L, 1);   % rand distribution
p = p/sum(p); % raw stochastic

[X, shifts] = generate_observations_Zn(x, N_arr(end), sigma, p); 

%% analytic moments

[M1, M2, M3] = Analytical_Moments_Zn(x, p);

%% test convergence
err1 = zeros(num_of_trials,1);
err2 = zeros(num_of_trials,1);
err3 = zeros(num_of_trials,1);
for j=1:num_of_trials
    N = N_arr(j);
    currentX = X(:,1:N);
    mu1 = mean(currentX,2);
    mu2 = (1/N)*currentX*(currentX')-sigma^2*eye(L);
   
    mu3 = zeros(L,L,L);
    for k=1:L
     %   for j=1:N
            mu3(:,:,k) = mu3(:,:,k) + currentX*diag(currentX(k,:))*currentX';
      %  end
    end
    mu3 = (1/N)*mu3;

   % A = kron(currentX,currentX)*currentX';
   
% ==>   mu3 = (1/N)*currentX*(currentX')-sigma^2*eye(L); % to ensure the bias

    err1(j) = norm(mu1(:)-M1(:));
    err2(j) = norm(mu2(:)-M2(:));
    err3(j) = norm(mu3(:)-M3(:));
end

%% plot
figure;
loglog(N_arr, err1);
grid on;
hold on;
loglog(N_arr,exp(log(err1(1))-.5*(log(N_arr)-log(N_arr(1)))),'r' );
legend('comvergence','expected rate')
title('First moment')

%%
figure
loglog(N_arr, err2);
hold on;
loglog(N_arr,exp(log(err2(1))-.5*(log(N_arr)-log(N_arr(1)))),'r' );
legend('comvergence','expected rate')
grid on;
title('Second moment')

%%
figure
loglog(N_arr, err3);
hold on;
loglog(N_arr,exp(log(err3(1))-.5*(log(N_arr)-log(N_arr(1)))),'r' );
legend('comvergence','expected rate')
grid on;
title('Third moment')

toc()