function [res] = B2(z)
    [m, n] = size(z);
    res = zeros(m, n);
    for i = 1:m
        for j = 1:n
            if (abs(z(i,j)) < 1e-8)
                res(i,j) = 1/6;
            else
                res(i,j) = (2 + z(i,j) + exp(z(i,j))*(-2 + z(i,j)))/((z(i,j))^3);
            end
        end
    end
end
