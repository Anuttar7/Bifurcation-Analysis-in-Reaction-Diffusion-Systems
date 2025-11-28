function [res] = B3(z)
    [m, n] = size(z);
    res = zeros(m, n);
    for i = 1:m
        for j = 1:n
            if (abs(z(i,j)) < 1e-8)
                res(i,j) = 1/6;
            else
                res(i,j) = (-4 - 3*z(i,j) - (z(i,j))^2 + exp(z(i,j))*(4 - z(i,j)))/((z(i,j))^3);
            end
        end
    end
end