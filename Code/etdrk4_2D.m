%Function that implements ETDRK4 method to solve Reaction-Diffusion System
%for 2D case
function [hist_u, hist_v] = etdrk4_2D(num_t, MU, PHI, PSI, u0, v0)
    %Parameters:
    D = 0.1;
    mu = MU;
    phi = PHI;
    psi = PSI;

    %Number of Iterations:
    Nt = num_t;

    %Size of Stencil:
    Nx = 50;
    Ny = 50;

    %delta variables:
    dx = 0.01;
    dy = 0.01;
    dt = 0.01;

    %Axes:
    x = dx:dx:Nx*dx;
    y = dy:dy:Ny*dy;
    t = dt:dt:Nt*dt;

    %Initialisation:
    u = zeros(Nx, Ny);
    v = zeros(Nx, Ny);
    
    for ii = 1:Nx
        for jj = 1:Ny
            u(ii,jj) = u0 - (2e-7)*(x(ii) - 0.1*y(jj) - 0.225)*(x(ii) - 0.1*y(jj) - 0.675);
            v(ii,jj) = v0 - (3e-5)*(x(ii) - 0.450) - (1.2e-4)*(y(jj) - 0.150);
        end
    end

    hist_u = zeros(Nx, Ny, Nt);
    hist_v = zeros(Nx, Ny, Nt);

    I = eye(2*Nx*Ny);

    %Creating Function Handles
    phi1f = @phi1;
    B1f = @B1;
    B2f = @B2;
    B3f = @B3;

    %Evaluating CF Parameters
    [zi1, ci1, r_inf1, temp, temp] = cf_parameters(phi1f, 12);
    [zi2, ci2, r_inf2, temp, temp] = cf_parameters(B1f, 12);
    [zi3, ci3, r_inf3, temp, temp] = cf_parameters(B2f, 12);
    [zi4, ci4, r_inf4, temp, temp] = cf_parameters(B3f, 12);

    %Calculating L and its exponentials
    L = calculate_L_2D(Nx, Ny, D, dx, dy);
    eL = expm(dt*L);
    eL2 = expm((dt/2)*L);
    
    %Main Loop
    for it = 1:Nt
        if (mod(it + 99, 5) == 0)
            fprintf('Iteration Number: %d\n', it);
        end

        %Transforming u and v into a single vector w
        w = reshape([reshape(u, Nx*Ny, 1), reshape(v, Nx*Ny, 1)], 2*Nx*Ny, 1);
        u = w(1:Nx*Ny);
        v = w(Nx*Ny + 1 : 2*Nx*Ny);

        %Evaluating ETDRK4 parameters
        a = eL2*w + dt*fn_approx_cf(zi1, ci1, r_inf1, N(u,v));
        a_u = a(1:Nx*Ny);
        a_v = a(Nx*Ny + 1 : 2*Nx*Ny);

        b = eL2*w + dt*fn_approx_cf(zi1, ci1, r_inf1, N(a_u, a_v));
        b_u = b(1:Nx*Ny);
        b_v = b(Nx*Ny + 1 : 2*Nx*Ny);
        
        c = eL2*a + dt*fn_approx_cf(zi1, ci1, r_inf1, 2*N(b_u, b_v) - N(u, v));
        c_u = c(1:Nx*Ny);
        c_v = c(Nx*Ny + 1 : 2*Nx*Ny);

        %Updating w
        w = eL*w + dt * (fn_approx_cf(zi2, ci2, r_inf2, N(u, v)) + fn_approx_cf(zi3, ci3, r_inf3, 2*(N(a_u, a_v) + N(b_u, b_v))) + fn_approx_cf(zi4, ci4, r_inf4, N(c_u, c_v)));

        %Extracting u and v from w
        u = reshape(real(w(1:Nx*Ny)), Ny, Nx);
        v = reshape(real(w(Nx*Ny + 1 : 2*Nx*Ny)), Ny, Nx);

        hist_u(:, :, it) = u;
        hist_v(:, :, it) = v;
    end

    %Function that returns Nonlinear part of the governing equation
    function [res] = N(u, v)
        res1 = u.*(1 - u) - (mu*(u.*v))./(u + phi);
        res2 = psi*v - (psi*(v.*v))./u;
        res = reshape([res1, res2], 2*Nx*Ny, 1);
    end

    %Function that approximates a given function using CF parameters
    function [yy] = fn_approx_cf(zi, ci, r_inf, F)
        yy = r_inf * F;
        for j = 1:12
            xx = (dt*L - zi(j)*I) \ F;
            yy = yy + ci(j) * xx;
        end
    end
end
