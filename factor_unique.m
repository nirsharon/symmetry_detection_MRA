% Helper function to get all positive divisors of n
function d = factor_unique(n)
    d = unique(factor_all(n));
end

function out = factor_all(n)
    out = [];
    for i = 1:n
        if mod(n, i) == 0
            out(end+1) = i;
        end
    end
end
