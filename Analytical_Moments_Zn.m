function [M1, M2, M3] = Analytical_Moments_Zn(x, rho)
% calculate the analytical moments according to formula
Cx = circulant(x);

M1 = Cx*rho;
M2 = Cx*diag(rho)*Cx';

n = length(x);

%--------------------------
direct_calc=0;
if direct_calc
    % M2 via a direct calculation
    m2 = zeros(n,n);
    for i=1:n
        for j=1:n
            for k=1:n
                m2(i,j) = m2(i,j) + x(1+mod(i-k,n))*x(1+mod(j-k,n))*rho(k);
            end
        end
    end

    % M3 via a direct calculation
    m3 = zeros(n,n,n);
    for i=1:n
        for j=1:n
            for l=1:n
                for k=1:n
                    m3(i,j,l) = m3(i,j,l) + x(1+mod(i-k,n))*x(1+mod(j-k,n))*x(1+mod(l-k,n))*rho(k);
                end
            end
        end
    end
    %
end
%--------------------------

M3 = zeros(n,n,n);
Crx= circulant(reverse(x));
for i=1:n
    M3(:,:,i) = (Cx*(diag(rho)*diag(Crx(:,i)))*Cx');
    % check:  norm(m3(:,:,i)-M3(:,:,i));
end



end