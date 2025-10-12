function divs = divisors_of_n(n)
    divs = [];
    for k = 1:n
        if mod(n, k) == 0
            divs(end+1) = k;
        end
    end
end