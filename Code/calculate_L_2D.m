function [L] = calculate_L_2D(Nx, Ny, D, dx, dy)
        L = zeros(2*Nx*Ny, 2*Nx*Ny);
        arrD = [1, D];

        ind = @(i,j,k) Nx*Ny*(k-1) + Ny*(i-1) + j;

        for i = 1:Nx
            for j = 1:Ny
                for k = 1:2
                    %Adding d2/dx2 terms
                    if (i == 1)
                        L(ind(i,j,k), ind(i+2,j,k)) = L(ind(i,j,k), ind(i+2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                        L(ind(i,j,k), ind(i+1,j,k)) = L(ind(i,j,k), ind(i+1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dx*dx));
                        L(ind(i,j,k), ind(i+1,j,k)) = L(ind(i,j,k), ind(i+1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i+2,j,k)) = L(ind(i,j,k), ind(i+2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                    elseif (i == 2)
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-1/(12*dx*dx));
                        L(ind(i,j,k), ind(i-1,j,k)) = L(ind(i,j,k), ind(i-1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dx*dx));
                        L(ind(i,j,k), ind(i+1,j,k)) = L(ind(i,j,k), ind(i+1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i+2,j,k)) = L(ind(i,j,k), ind(i+2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                    elseif (i == Nx - 1)
                        L(ind(i,j,k), ind(i-2,j,k)) = L(ind(i,j,k), ind(i-2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                        L(ind(i,j,k), ind(i-1,j,k)) = L(ind(i,j,k), ind(i-1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dx*dx));
                        L(ind(i,j,k), ind(i+1,j,k)) = L(ind(i,j,k), ind(i+1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-1/(12*dx*dx));
                    elseif (i == Nx)
                        L(ind(i,j,k), ind(i-2,j,k)) = L(ind(i,j,k), ind(i-2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                        L(ind(i,j,k), ind(i-1,j,k)) = L(ind(i,j,k), ind(i-1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dx*dx));
                        L(ind(i,j,k), ind(i-1,j,k)) = L(ind(i,j,k), ind(i-1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i-2,j,k)) = L(ind(i,j,k), ind(i-2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                    else
                        L(ind(i,j,k), ind(i-2,j,k)) = L(ind(i,j,k), ind(i-2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                        L(ind(i,j,k), ind(i-1,j,k)) = L(ind(i,j,k), ind(i-1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dx*dx));
                        L(ind(i,j,k), ind(i+1,j,k)) = L(ind(i,j,k), ind(i+1,j,k)) + arrD(k)*(16/(12*dx*dx));
                        L(ind(i,j,k), ind(i+2,j,k)) = L(ind(i,j,k), ind(i+2,j,k)) + arrD(k)*(-1/(12*dx*dx));
                    end

                    %Adding d2/dy2 terms
                    if (j == 1)
                        L(ind(i,j,k), ind(i,j+2,k)) = L(ind(i,j,k), ind(i,j+2,k)) + arrD(k)*(-1/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+1,k)) = L(ind(i,j,k), ind(i,j+1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+1,k)) = L(ind(i,j,k), ind(i,j+1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+2,k)) = L(ind(i,j,k), ind(i,j+2,k)) + arrD(k)*(-1/(12*dy*dy));
                    elseif (j == 2)
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-1/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j-1,k)) = L(ind(i,j,k), ind(i,j-1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+1,k)) = L(ind(i,j,k), ind(i,j+1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+2,k)) = L(ind(i,j,k), ind(i,j+2,k)) + arrD(k)*(-1/(12*dy*dy));
                    elseif (j == Ny - 1)
                        L(ind(i,j,k), ind(i,j-2,k)) = L(ind(i,j,k), ind(i,j-2,k)) + arrD(k)*(-1/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j-1,k)) = L(ind(i,j,k), ind(i,j-1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+1,k)) = L(ind(i,j,k), ind(i,j+1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-1/(12*dy*dy));
                    elseif (j == Ny)
                        L(ind(i,j,k), ind(i,j-2,k)) = L(ind(i,j,k), ind(i,j-2,k)) + arrD(k)*(-1/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j-1,k)) = L(ind(i,j,k), ind(i,j-1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j-1,k)) = L(ind(i,j,k), ind(i,j-1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j-2,k)) = L(ind(i,j,k), ind(i,j-2,k)) + arrD(k)*(-1/(12*dy*dy));
                    else
                        L(ind(i,j,k), ind(i,j-2,k)) = L(ind(i,j,k), ind(i,j-2,k)) + arrD(k)*(-1/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j-1,k)) = L(ind(i,j,k), ind(i,j-1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j,k))   = L(ind(i,j,k), ind(i,j,k))   + arrD(k)*(-30/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+1,k)) = L(ind(i,j,k), ind(i,j+1,k)) + arrD(k)*(16/(12*dy*dy));
                        L(ind(i,j,k), ind(i,j+2,k)) = L(ind(i,j,k), ind(i,j+2,k)) + arrD(k)*(-1/(12*dy*dy));
                    end
                end
            end
        end
    end
