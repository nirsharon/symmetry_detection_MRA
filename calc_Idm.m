function [Idm] = calc_Idm(d, m)
% This function outputs all the expected entries of the power spectrum that
% are zerod-out in symmetry Zm

% OLD VERSION
%Idm = zeros(d,1);
%for j=0:(d-1)
%    if mod(j,m)==0
%        Idm(j+1)=1;
%    end
%end

% NEW VERSION: 
% 1 if multiple of m, 0 otherwise
Idm = (mod(0:(d-1), m) == 0)';

end