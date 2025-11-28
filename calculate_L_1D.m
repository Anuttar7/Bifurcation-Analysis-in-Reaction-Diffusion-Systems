function [L] = calculate_L_1D(Nx, D, dx)
        L = zeros(2*Nx, 2*Nx);
        arrD = [1, D];

        ind = @(i,k) Nx*(k-1) + i;

        for i = 1:Nx
            for k = 1:2
                %Adding d2/dx2 terms
                if (i == 1)
                    L(ind(i,k), ind(i+2,k)) = L(ind(i,k), ind(i+2,k)) + arrD(k)*(-1/(12*dx*dx));
                    L(ind(i,k), ind(i+1,k)) = L(ind(i,k), ind(i+1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-30/(12*dx*dx));
                    L(ind(i,k), ind(i+1,k)) = L(ind(i,k), ind(i+1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i+2,k)) = L(ind(i,k), ind(i+2,k)) + arrD(k)*(-1/(12*dx*dx));
                elseif (i == 2)
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-1/(12*dx*dx));
                    L(ind(i,k), ind(i-1,k)) = L(ind(i,k), ind(i-1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-30/(12*dx*dx));
                    L(ind(i,k), ind(i+1,k)) = L(ind(i,k), ind(i+1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i+2,k)) = L(ind(i,k), ind(i+2,k)) + arrD(k)*(-1/(12*dx*dx));
                elseif (i == Nx - 1)
                    L(ind(i,k), ind(i-2,k)) = L(ind(i,k), ind(i-2,k)) + arrD(k)*(-1/(12*dx*dx));
                    L(ind(i,k), ind(i-1,k)) = L(ind(i,k), ind(i-1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-30/(12*dx*dx));
                    L(ind(i,k), ind(i+1,k)) = L(ind(i,k), ind(i+1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-1/(12*dx*dx));
                elseif (i == Nx)
                    L(ind(i,k), ind(i-2,k)) = L(ind(i,k), ind(i-2,k)) + arrD(k)*(-1/(12*dx*dx));
                    L(ind(i,k), ind(i-1,k)) = L(ind(i,k), ind(i-1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-30/(12*dx*dx));
                    L(ind(i,k), ind(i-1,k)) = L(ind(i,k), ind(i-1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i-2,k)) = L(ind(i,k), ind(i-2,k)) + arrD(k)*(-1/(12*dx*dx));
                else
                    L(ind(i,k), ind(i-2,k)) = L(ind(i,k), ind(i-2,k)) + arrD(k)*(-1/(12*dx*dx));
                    L(ind(i,k), ind(i-1,k)) = L(ind(i,k), ind(i-1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i,k))   = L(ind(i,k), ind(i,k))   + arrD(k)*(-30/(12*dx*dx));
                    L(ind(i,k), ind(i+1,k)) = L(ind(i,k), ind(i+1,k)) + arrD(k)*(16/(12*dx*dx));
                    L(ind(i,k), ind(i+2,k)) = L(ind(i,k), ind(i+2,k)) + arrD(k)*(-1/(12*dx*dx));
                end
            end
        end
    end