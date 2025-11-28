function [res] = phi1(z)
    [m, n] = size(z);
    res = zeros(m, n);
    for i = 1:m
        for j = 1:n
            if (abs(z(i,j)) < 1e-6)
                res(i,j) = 1;
            else
                res(i,j) = (exp(z(i,j)) - 1)/z(i,j);
            end
        end
    end
end
