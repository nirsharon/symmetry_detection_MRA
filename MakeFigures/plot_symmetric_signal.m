% script name: "plot_symmetric_signal"

% cleaning the table
close all;
clear;

% basic parameters
n = 60;
m = 12;

to_save = 1;

% the signal
%f = @(x) sin(140 * pi * x);
%base_x = linspace(-1.5, 1.5, n);
base_x = randn(n,1);
x = get_symmetric_Zm_signal_Zn(n, n/m, base_x);

% the distribution
rho = abs(rand(n,1));
rho = rho/sum(rho);

% the analytical moments (just for comparison)
Cx = circulant(x);
M1 = Cx*rho;
M2 = Cx*diag(rho)*Cx';
F = fft(eye(n));

if to_save
    folder_name = 'sym_signal';
    if isfolder(folder_name)
        cd(folder_name);
    else
        mkdir(folder_name)
        cd(folder_name)
    end

    figure;
    plot(x,'LineWidth',1.8);

    name_it = ['d_', num2str(n),'m_', num2str(m),'_the_signal' ] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    print('-depsc2',name_it);
    print('-depsc2',name_it);

    figure;
    stem(real(diag(F*M2*F')), 'filled','LineWidth',1.4)
    title('The power spectrum')

    name_it = ['d_', num2str(n),'m_', num2str(m),'_the_PS' ] ;
    saveas(gcf, name_it ,'fig');
    saveas(gcf,name_it,'jpg');
    print('-depsc2',name_it);
    print('-depsc2',name_it);

    cd '../'
else
    figure;
    plot(x,'LineWidth',1.5);
    figure;
    stem(real(diag(F*M2*F')), 'filled','LineWidth',1.2)
    title('The power spectrum')
end
