function [phi_hat, eta_hat] = evolve_wave(phi_hat, eta_hat,t_in,p)
% Evolves wave field in F-space through nsteps of size dt

t = t_in;

%{
for n=1:p.nsteps_impact  
    
    [rhs1_1, rhs2_1] = compute_rhs_full(phi_hat, eta_hat, t,p); 
    [rhs1_2, rhs2_2] = compute_rhs_full(phi_hat+p.dt/2.*rhs1_1, eta_hat+p.dt/2.*rhs2_1,t+p.dt/2,p); 
    [rhs1_3, rhs2_3] = compute_rhs_full(phi_hat+p.dt/2.*rhs1_2, eta_hat+p.dt/2.*rhs2_2,t+p.dt/2,p); 
    [rhs1_4, rhs2_4] = compute_rhs_full(phi_hat+p.dt.*rhs1_3, eta_hat+p.dt.*rhs2_3,t+p.dt,p);
    phi_hat = phi_hat + p.dt/6*(rhs1_1 + 2*rhs1_2 + 2*rhs1_3 + rhs1_4);
    eta_hat = eta_hat + p.dt/6*(rhs2_1 + 2*rhs2_2 + 2*rhs2_3 + rhs2_4);

    t = t+p.dt;
end      
%}

eps = p.dt/20;
phi_hat0 = phi_hat;
eta_hat0 = eta_hat;
phi_hat1 = 0;
phi_hat2 = 0;
phi_hat3 = 0;
eta_hat1 = 0;
eta_hat2 = 0;
eta_hat3 = 0;

for n=1:min(30,10*p.nsteps_impact)
    
    RHS1 = - p.g(t)*eta_hat + p.Bo*p.Ky2_new*eta_hat + p.Bo*eta_hat*p.Kx2_new.';
    phi_hat = (speye(p.Ny,p.Ny) - eps*2*p.nu0*p.Ky2_new)\(phi_hat*(speye(p.Nx,p.Nx) + eps*2*p.nu0*p.Kx2_new.') + eps*RHS1);
    phi_hat = ((speye(p.Ny,p.Ny) + eps*2*p.nu0*p.Ky2_new.')*phi_hat + eps*RHS1)/(speye(p.Nx,p.Nx) - eps*2*p.nu0*p.Kx2_new.');

    RHS2 = DtN(phi_hat,p);
    eta_hat = (speye(p.Ny,p.Ny) - eps*2*p.nu0*p.Ky2_new)\(eta_hat*(speye(p.Nx,p.Nx) + eps*2*p.nu0*p.Kx2_new.') + eps*RHS2);
    eta_hat = ((speye(p.Ny,p.Ny) + eps*2*p.nu0*p.Ky2_new)*eta_hat + eps*RHS2)/(speye(p.Nx,p.Nx) - eps*2*p.nu0*p.Kx2_new.');

    t = t+p.dt;
    
    if ( n == 10 )
        phi_hat1 = phi_hat;
        eta_hat1 = eta_hat;
    elseif n == 20
        phi_hat2 = phi_hat;
        eta_hat2 = eta_hat;
    elseif n == 30
        phi_hat3 = phi_hat;
        eta_hat3 = eta_hat;
    end
end

eps = p.dt/2;
%surf(p.xx,p.yy,real(ifft2(phi_hat)))

for n=4:p.nsteps_impact  
    
    RHS1 = - p.g(t)*eta_hat + p.Bo*p.Ky2_new*eta_hat + p.Bo*eta_hat*p.Kx2_new.';
    phi_hat = (speye(p.Ny,p.Ny) - eps*2*p.nu0*p.Ky2_new)\(phi_hat*(speye(p.Nx,p.Nx) + eps*2*p.nu0*p.Kx2_new.') + eps*RHS1);
    phi_hat = ((speye(p.Ny,p.Ny) + eps*2*p.nu0*p.Ky2_new.')*phi_hat + eps*RHS1)/(speye(p.Nx,p.Nx) - eps*2*p.nu0*p.Kx2_new.');
    
    
    

    RHS2 = DtN(phi_hat,p);
    eta_hat = (speye(p.Ny,p.Ny) - eps*2*p.nu0*p.Ky2_new)\(eta_hat*(speye(p.Nx,p.Nx) + eps*2*p.nu0*p.Kx2_new.') + eps*RHS2);
    eta_hat = ((speye(p.Ny,p.Ny) + eps*2*p.nu0*p.Ky2_new)*eta_hat + eps*RHS2)/(speye(p.Nx,p.Nx) - eps*2*p.nu0*p.Kx2_new.');

    t = t+p.dt;
    
    clf;
    %surf(p.xx,p.yy,abs(phi_hat))
    %hold on;
    surf(p.xx,p.yy,real(ifft2(eta_hat)))
    %axis([-24 24 -24 24 -1 1])
    colormap summer
    shading interp
    pause(0.2)
end      
%rhs1      =   -p.g(t).*eta_hat + p.Bo*p.K2.*eta_hat + 2*p.nu0.*p.K2.*phi_hat;
%rhs2 = DtN(phi_hat,p) + 2*p.nu0.*p.K2.*eta_hat; % DtN computes phi_z_hat    
end

     
