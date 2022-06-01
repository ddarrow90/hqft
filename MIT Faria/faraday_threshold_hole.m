close all
clc
clear all

add_to_path

% The threshold test should be run with the SAME Nx and Lx as the
% simulations of interest. Also make sure that the geometry is the same
Gam = [0.990:0.002:1.010]*3.82; 
Nx = 512; Ny = Nx; 
Lx = 48; Ly = Lx; dt_desired = min(Lx/Nx,Ly/Ny)/10;

parfor n=1:length(Gam)
    p = problem_setup_hole(Nx,Ny,Lx,Ly,Gam(n),dt_desired);

    %Define bottom
    p.Dhole = 13e-3/p.xF; %13mm = 2.18*\lambda_F
    
    deep = (sqrt(p.yy.^2 + p.xx.^2)< (p.Dhole./2));

    p.d = p.d0_deep.*deep + p.d0_shallow.*(~deep);
    p.h = p.h0_deep.*deep + p.h0_shallow.*(~deep);
    
    % Number of impacts (time to wait to check the threshold)
    p.nimpacts = 500; %1000 is a good number 
    [is_stable(n),eta_max(n)] = faraday_threshold(p); 
end

% save('stability_diagram.mat','a','is_stable','eta_max(n)');
plot(Gam,eta_max,'rx'); xlabel('Gamma'); ylabel('max(eta)');

saveas(gcf,['Faraday_threshold_N_',num2str(Nx),'.fig']);




