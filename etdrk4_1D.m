%Function that implements ETDRK4 method to solve Reaction-Diffusion System
%for 1D case
function [hist_u, hist_v] = etdrk4_1D(num_t, MU, PHI, PSI, u0, v0)
    %Parameters:
    D = 0.1;
    mu = MU;
    phi = PHI;
    psi = PSI;

    %Number of Iterations:
    Nt = num_t;

    %Size of Stencil:
    Nx = 50;

    %delta variables:
    dx = 0.1;
    dt = 0.1;

    %Axes:
    x = dx:dx:Nx*dx;
    t = dt:dt:Nt*dt;

    %Initialisation:
    u = zeros(Nx, 1);
    v = zeros(Nx, 1);
    
    for ii = 1:Nx
        u(ii) = u0 - (2e-7)*(x(ii) - 2.25)*(x(ii) - 6.75);
        v(ii) = v0 - (3e-5)*(x(ii) - 4.50);
    end

    hist_u = zeros(Nx, Nt);
    hist_v = zeros(Nx, Nt);

    %Defining Identity Matrix
    I = eye(2*Nx);

    %Defining Function Handles
    phi1f = @phi1;
    B1f = @B1;
    B2f = @B2;
    B3f = @B3;

    %Obtaining CF Parameters
    [zi1, ci1, r_inf1, temp, temp] = cf_parameters(phi1f, 12);
    [zi2, ci2, r_inf2, temp, temp] = cf_parameters(B1f, 12);
    [zi3, ci3, r_inf3, temp, temp] = cf_parameters(B2f, 12);
    [zi4, ci4, r_inf4, temp, temp] = cf_parameters(B3f, 12);

    %Calculating L and its exponetials
    L = calculate_L_1D(Nx, D, dx);
    eL = expm(dt*L);
    eL2 = expm((dt/2)*L);
    
    %Main Loop
    for it = 1:Nt
        if (mod(it + 99, 100) == 0)
            fprintf('Iteration Number: %d\n', it);
        end

        %Reshaping all variables in a single vector
        w = reshape([u, v], 2*Nx, 1);
        u = w(1:Nx);
        v = w(Nx + 1 : end);

        %Evaluating ETDRK4 Parameters
        a = eL2*w + dt*fn_approx_cf(zi1, ci1, r_inf1, N(u,v));
        a_u = a(1:Nx);
        a_v = a(Nx + 1 : end);

        b = eL2*w + dt*fn_approx_cf(zi1, ci1, r_inf1, N(a_u, a_v));
        b_u = b(1:Nx);
        b_v = b(Nx + 1 : 2*Nx);
        
        c = eL2*a + dt*fn_approx_cf(zi1, ci1, r_inf1, 2*N(b_u, b_v) - N(u, v));
        c_u = c(1:Nx);
        c_v = c(Nx + 1 : 2*Nx);

        %Updating w
        w = eL*w + dt * (fn_approx_cf(zi2, ci2, r_inf2, N(u, v)) + fn_approx_cf(zi3, ci3, r_inf3, 2*(N(a_u, a_v) + N(b_u, b_v))) + fn_approx_cf(zi4, ci4, r_inf4, N(c_u, c_v)));

        %Extracting u and v from w
        u = real(w(1:Nx));
        v = real(w(Nx + 1 : end));

        hist_u(:, it) = u;
        hist_v(:, it) = v;
    end

    %Function that returns non-linear part of equation
    function [res] = N(u, v)
        res1 = u.*(1 - u) - (mu*(u.*v))./(u + phi);
        res2 = psi*v - (psi*(v.*v))./u;
        res = reshape([res1, res2], 2*Nx, 1);
    end 

    %Function that returns cf approximation to a given function
    function [yy] = fn_approx_cf(zi, ci, r_inf, F)
        yy = r_inf * F;
        for j = 1:12
            xx = (dt*L - zi(j)*I) \ F;
            yy = yy + ci(j) * xx;
        end
    end
end