function [phi_hat, eta_hat] = evolve_wave_BDF4(phi_hat, eta_hat,t_in,p, type)
% Evolves wave field in F-space through nsteps of size dt
debug0 = 0;

if type == 0
    [phi_hat,eta_hat] = evolve_wave(phi_hat,eta_hat,t_in,p);
    return
end

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
phi_hat0 = phi_hat;
eta_hat0 = eta_hat;
phi_hat1 = 0;
phi_hat2 = 0;
eta_hat1 = 0;
eta_hat2 = 0;

temp_vec = phi_hat(:); %[vec(phi_hat); vec(eta_hat)];

for n=1:min(2*smallN,smallN*p.nsteps_impact)
    
    RHS1 = - p.g(t)*eta_hat + p.Bo*p.Ky2_new*eta_hat + p.Bo*eta_hat*p.Kx2_new.';
    for x=1:p.Nx
        phi_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)\(phi_hat(:,x) + eps*RHS1(:,x));
    end

    RHS2 = DtN_new(phi_hat,p);
    for x=1:p.Nx
        eta_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)\(eta_hat(:,x) + eps*RHS2(:,x));
    end
    t = t+p.dt;
    
    if ( n == smallN && type > 1)
        phi_hat1 = phi_hat;
        eta_hat1 = eta_hat;
        plot_sol(eta_hat1,p,0, 1*p.dt)
    elseif n == 2*smallN && type > 1
        phi_hat2 = phi_hat;
        eta_hat2 = eta_hat;
        plot_sol(eta_hat2,p,0, 2*p.dt)
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

L = kron((speye(p.Nx,p.Nx) - eps*2*p.nu0*p.Kx2_new),speye(p.Ny,p.Ny))...
    + kron(speye(p.Nx,p.Nx),- eps*2*p.nu0*p.Ky2_new);
I_big = kron(speye(p.Nx,p.Nx),speye(p.Ny,p.Ny));
D2_big = kron(p.Kx2_new,speye(p.Ny,p.Ny)) + kron(speye(p.Nx,p.Nx),p.Ky2_new);
L_big = I_big - eps*2*p.nu0*D2_big;

O_big = kron(sparse(p.Nx,p.Nx),sparse(p.Ny,p.Ny));

M_big = [L_big, -p.Bo*D2_big; O_big, L_big];

for n=3:p.nsteps_impact  
    
    RHS1 = - p.g(t)*eta_hat + p.Bo*p.Ky2_new*eta_hat + p.Bo*eta_hat*p.Kx2_new.';
    if type > 1
        RHS2 = c1*phi_hat + c2*phi_hat2 + c3*phi_hat1 + c4*phi_hat0;
        phi_hat0 = phi_hat1;
        phi_hat1 = phi_hat2;
        phi_hat2 = phi_hat;
    else
        RHS2 = phi_hat;
    end
    temp_vec = L\(eps*RHS1(:) + RHS2(:));
    phi_hat = reshape(temp_vec, p.Ny, p.Nx);
    %for x=1:p.Nx
        
    %    phi_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)...
    %                   \(eps*RHS1(:,x) + RHS2(:,x));
    %end
    
    RHS1 = DtN_new(phi_hat,p);
    if type > 1
        RHS2 = c1*eta_hat + c2*eta_hat2 + c3*eta_hat1 + c4*eta_hat0;
        eta_hat0 = eta_hat1;
        eta_hat1 = eta_hat2;
        eta_hat2 = eta_hat;
    else
        RHS2 = eta_hat;
    end
    
    temp_vec = L\(eps*RHS1(:) + RHS2(:));
    eta_hat = reshape(temp_vec, p.Ny, p.Nx);

    %for x=1:p.Nx
    %    eta_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)...
    %                    \(eps*RHS1(:,x) + RHS2(:,x));
    %end
    t = t+p.dt;
    plot_sol(eta_hat,p,0, n*p.dt)
end      
%rhs1      =   -p.g(t).*eta_hat + p.Bo*p.K2.*eta_hat + 2*p.nu0.*p.K2.*phi_hat;
%rhs2 = DtN(phi_hat,p) + 2*p.nu0.*p.K2.*eta_hat; % DtN computes phi_z_hat    
end

     
