function x_out = reverse(x)
% x might be a matrix, and then each column is reversed 
x_out = zeros(size(x));
x_out(1,:) = x(1,:);
x_out(2:end,:) = flipud(x(2:end,:));
end