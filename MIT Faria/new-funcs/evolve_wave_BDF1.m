function [phi_hat, eta_hat] = evolve_wave_BDF1(phi_hat, eta_hat,t_in,p)
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
eps = p.dt;
%surf(p.xx,p.yy,real(ifft2(phi_hat)))
for n=1:p.nsteps_impact  
    
    RHS1 = - p.g(t)*eta_hat + p.Bo*p.Ky2_new*eta_hat + p.Bo*eta_hat*p.Kx2_new.';
    for x=1:p.Nx
        phi_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)\(phi_hat(:,x) + eps*RHS1(:,x));
    end
    
    

    RHS2 = DtN(phi_hat,p);
    for x=1:p.Nx
        eta_hat(:,x) = (speye(p.Ny,p.Ny)*(1 - eps*2*p.nu0*p.Kx2_new(x,x)) - eps*2*p.nu0*p.Ky2_new)\(eta_hat(:,x) + eps*RHS2(:,x));
    end

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

     
