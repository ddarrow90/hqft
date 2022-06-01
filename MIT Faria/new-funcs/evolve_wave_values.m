function [phi_hat, eta_hat] = evolve_wave_values(phi_hat, eta_hat,t_in,p, type)
% Evolves wave field in F-space through nsteps of size dt
debug0 = 1;

if type == 0
    [phi_hat,eta_hat] = evolve_wave(phi,eta_hat,t_in,p);
    return
end

phi = fft2(phi_hat);
eta = fft2(eta_hat);

t = t_in;

    function plot_sol(ifield,p,fixaxes, steps)
        if debug0 == 0
            return
        end
        clf;
        %surf(p.xx,p.yy,abs(phi_hat))
        %hold on;
        
        %axis([-24 24 -24 24 -1 1])
        surf(p.xx,p.yy,real(ifft2(ifield)))
        if fixaxes == 1
            axis([-24 24 -24 24 -1e-3 1e-3])
        end
        colormap summer
        shading interp
        title(["time = ", num2str(steps)])
        pause(0.2)
    end
smallN = 10;
eps = p.dt/smallN;
phi_0 = phi;
eta_0 = eta;
phi_1 = 0;
phi_2 = 0;
eta_1 = 0;
eta_2 = 0;

temp_vec = [phi(:); eta(:)]; %[vec(phi_hat); vec(eta_hat)];
L = kron((speye(p.Nx,p.Nx) - eps*2*p.nu0*p.D2x_new),speye(p.Ny,p.Ny))...
    + kron(speye(p.Nx,p.Nx),- eps*2*p.nu0*p.D2y_new);
I_big = kron(speye(p.Nx,p.Nx),speye(p.Ny,p.Ny));
D2_big = kron(p.D2x_new,speye(p.Ny,p.Ny)) + kron(speye(p.Nx,p.Nx),p.D2y_new);
L_big = I_big - eps*2*p.nu0*D2_big;

O_big = kron(sparse(p.Nx,p.Nx),sparse(p.Ny,p.Ny));

M_big = [L_big, -p.Bo*eps*D2_big; O_big, L_big];

for n=1:min(2*smallN,smallN*p.nsteps_impact)
    
    RHS1_1 = - 0*p.g(t)*eta;
    RHS2_1 = phi;

    RHS1_2 = ifft2(DtN_new(fft2(phi),p));
    RHS2_2 = eta;
    

    Mtime_big = [O_big, eps*p.g(t)*I_big; O_big, O_big];

    temp_vec = (M_big + Mtime_big)\(eps*[RHS1_1(:); RHS1_2(:)] + [RHS2_1(:); RHS2_2(:)]);
    phi = reshape(temp_vec(1:end/2).', p.Ny, p.Nx);
    eta = reshape(temp_vec(end/2+1:end).', p.Ny, p.Nx);

    t = t+p.dt/smallN;
    if ( n == smallN && type > 1)
        phi_1 = phi;
        eta_1 = eta;
        plot_sol(eta_1,p,0, 1*p.dt)
    elseif n == 2*smallN && type > 1
        phi_2 = phi;
        eta_2 = eta;
        plot_sol(eta_2,p,0, 2*p.dt)
    end
end

if type == 4
    eps = (12/25)*p.dt;
    c1 = 48/25;
    c2 = -36/25;
    c3 = 16/25;
    c4 = -3/25;
elseif type == 3
    eps = (6/11)*p.dt;
    c1 = 18/11;
    c2 = -9/11;
    c3 = 2/11;
    c4 = 0;
elseif type == 2
    eps = (2/3)*p.dt;
    c1 = 4/3;
    c2 = -1/3;
    c3 = 0;
    c4 = 0;
else
    eps = p.dt;
    c1 = 1;
    c2 = 0;
    c3 = 0;
    c4 = 0;
end



for n=3:p.nsteps_impact  
    
    RHS1_1 = - 0*p.g(t)*eta;
    if type > 1
        RHS2_1 = c1*phi + c2*phi_2 + c3*phi_1 + c4*phi_hat0;
        phi_hat0 = phi_1;
        phi_1 = phi_2;
        phi_2 = phi;
    else
        RHS2_1 = phi;
    end


    RHS1_2 = ifft2(DtN_new(fft2(phi),p));
    if type > 1
        RHS2_2 = c1*eta + c2*eta_2 + c3*eta_1 + c4*eta_0;
        eta_0 = eta_1;
        eta_1 = eta_2;
        eta_2 = eta;
    else
        RHS2_2 = eta;
    end
    Mtime_big = [O_big, eps*p.g(t)*I_big; O_big, O_big];

    temp_vec = (M_big + Mtime_big)\(eps*[RHS1_1(:); RHS1_2(:)] + [RHS2_1(:); RHS2_2(:)]);
    phi = reshape(temp_vec(1:end/2).', p.Ny, p.Nx);
    eta = reshape(temp_vec(end/2+1:end).', p.Ny, p.Nx);

    
    %for x=1:p.Nx
        
    %    phi_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)...
    %                   \(eps*RHS1(:,x) + RHS2(:,x));
    %end
    

    %for x=1:p.Nx
    %    eta_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)...
    %                    \(eps*RHS1(:,x) + RHS2(:,x));
    %end
    t = t+p.dt;
    plot_sol(eta,p,0, n*p.dt)
end      
%rhs1      =   -p.g(t).*eta_hat + p.Bo*p.K2.*eta_hat + 2*p.nu0.*p.K2.*phi_hat;
%rhs2 = DtN(phi_hat,p) + 2*p.nu0.*p.K2.*eta_hat; % DtN computes phi_z_hat    
if debug0 == 1
    pause()
end
eta_hat = fft2(eta);
phi_hat = fft2(phi);
end

     
