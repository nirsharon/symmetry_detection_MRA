function d = find_divisors(n)
    d = find(~mod(n, 1:n));
end