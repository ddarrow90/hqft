% Delete unnecessary variables to reduce result file size

%% Load files
cd('results')
files = dir('*.mat');

for file = files'

        file.name
        %load(file.name,'p','x_data','y_data','t_data','theta','wavefield');
        load(file.name,'p','x_data','y_data','t_data','theta');
        
        
        p2.Bo = p.Bo;
        p2.DtN_method = p.DtN_method;
        p2.Gam = p.Gam;
        p2.Gamma_max_deep = p.Gamma_max_deep;
        p2.Gamma_max_shallow = p.Gamma_max_shallow;
        p2.Lx = p.Lx;
        p2.Ly = p.Ly;
        p2.M = p.M;
        p2.Nx = p.Nx;
        p2.Ny = p.Ny;
        p2.Reynolds = p.Reynolds;
        p2.TF = p.TF;
        p2.c4 = p.c4;
        p2.cf_air = p.cf_air;
        p2.cf_impact = p.cf_impact;
        p2.d0_deep = p.d0_deep;
        p2.d0_shallow = p.d0_shallow;
        p2.d_deep = p.d_deep;
        p2.d_shallow = p.d_shallow;
        p2.drop_density = p.drop_density;
        p2.drop_mass = p.drop_mass;
        p2.drop_radius = p.drop_radius;
        p2.drop_type = p.drop_type;
        p2.dt = p.dt;
        p2.dt_desired = p.dt_desired;
        p2.g0 = p.g0;
        p2.h0_deep = p.h0_deep;
        p2.h0_shallow = p.h0_shallow;
        p2.lambdaf_deep = p.lambdaf_deep;
        p2.lambdaf_shallow = p.lambdaf_shallow;
        p2.mu_air = p.mu_air;
        p2.nsteps_impact = p.nsteps_impact;
        p2.nu = p.nu;
        p2.nu0= p.nu0;
        p2.omega0 = p.omega0;
        p2.rho = p.rho;
        p2.sig = p.sig;
        p2.store_wavefield = p.store_wavefield;
        p2.w0 = p.w0;
        p2.xF = p.xF;
        p2.xx = p.xx;
        p2.yy = p.yy;
        p2.xi = p.xi;
        p2.yi = p.yi;
        p2.ui = p.ui;
        p2.vi = p.vi;
        p2.nimpacts = p.nimpacts;
        p2.Dhole = p.Dhole;
%        p2.L = p.L;
        p2.d = p.d;
        
        clearvars p
        
        p = p2;
        
    
        % Delete old file
        delete(file.name)
        % Save lighter variables
        %save(file.name,'p','x_data','y_data','t_data','theta','wavefield');
         save(file.name,'p','x_data','y_data','t_data','theta');
                
end


cd('../')
