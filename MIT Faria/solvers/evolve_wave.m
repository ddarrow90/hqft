function [phi_hat, eta_hat,psi_hat] = evolve_wave(phi_hat, eta_hat,psi_hat,t_in,p)
% Evolves wave field in F-space through nsteps of size dt

t = t_in;

for n=1:p.nsteps_impact  
    
    [rhs1_1, rhs2_1, rhs3_1] = compute_rhs_full(phi_hat, eta_hat, psi_hat, t,p); 
    [rhs1_2, rhs2_2, rhs3_2] = compute_rhs_full(phi_hat+p.dt/2.*rhs1_1, eta_hat+p.dt/2.*rhs2_1,psi_hat+p.dt/2.*rhs3_1, t+p.dt/2,p); 
    [rhs1_3, rhs2_3, rhs3_3] = compute_rhs_full(phi_hat+p.dt/2.*rhs1_2, eta_hat+p.dt/2.*rhs2_2,psi_hat+p.dt/2.*rhs3_2, t+p.dt/2,p); 
    [rhs1_4, rhs2_4, rhs3_4] = compute_rhs_full(phi_hat+p.dt.*rhs1_3, eta_hat+p.dt.*rhs2_3,psi_hat+p.dt.*rhs3_3, t+p.dt,p);
    phi_hat = phi_hat + p.dt/6*(rhs1_1 + 2*rhs1_2 + 2*rhs1_3 + rhs1_4);
    eta_hat = eta_hat + p.dt/6*(rhs2_1 + 2*rhs2_2 + 2*rhs2_3 + rhs2_4);
    %psi_hat = eta_hat + p.dt/6*(rhs3_1 + 2*rhs3_2 + 2*rhs3_3 + rhs3_4);

    t = t+p.dt;
end      

    
end

     
