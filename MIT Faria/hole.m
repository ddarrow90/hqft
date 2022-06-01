function hole(xi,yi,ui,vi,no_bottom, mem, system)
%% questions i have

% Wavenumbers are currently stored [ 1:N, 0, -N:-1]. wtf?
% DONE


% matrix multiplication is currently stored in super inefficient way with
% dense matrix of wavenumbers
% DONE


% we spend a long time computing unnecessary IFFTs in compute_slope_eta:
%  - we could use a point-cloud method
%  - we could just IFFT the specific point
% DONE


% note: the nearest point code can be redone using linear evaluation
% operators
% DONE


% in drop_impact, we add a delta function to simulate a particle impact;
% our spectral approach assumes smoothness, though, and we also get drastic
% aliasing issues when the particle is between lattice points. Consider
% smoothing between a few points


% are we sure we want such a singular contribution with each drop? at the
% very least, it's being resolved very poorly


% the wave equation is linear, guys, let's not use RK4
% DONE


% we use a first-order method in time for the droplets, which kinda ruins
% the whole RK4 thing. Since the velocity is a linear function of the known
% field, we can use an implicit method here as well. Or just RK4 here?

%% example setup
%xi = -22.3;
%yi = 0;
%ui = 0.7;
%vi = 0.2;
convert_ICs_to_lshallow = 0;
%no_bottom = 0;
spin = 0;
% Set parameters

Gam = mem*3.815;
Nx = 512; Ny = Nx; % 512
Lx = 48; Ly = Lx; dt_desired = min(Lx/Nx,Ly/Ny)/10;

p = problem_setup_hole(Nx,Ny,Lx,Ly,Gam,dt_desired);

% Drop's initial position

% Refer IP with respect to \lamba_shallow
if convert_ICs_to_lshallow
    xi = xi.*p.lambdaf_shallow./p.lambdaf_deep;
    yi = yi.*p.lambdaf_shallow./p.lambdaf_deep;
end

p.xi = xi; p.yi = yi; p.ui= ui; p.vi = vi; p.spin = spin;
p.nimpacts = 5001;      % Number of impacts
%p.dt = 0.008;
p.nsteps_impact = 107;round(107*0.01/p.dt); % 110, 160


width = 0.001;
type = 0;
plotoption = 10;

if type > 0
%    p.cf_air = p.cf_air*0.1;
%    p.G = 5*p.G;
%    p.cf_impact = p.cf_impact*4;

%    norm_h = 0.25;
 %   p.d0_deep = 15*p.d0_deep; % below 0.008, above 0.0002
 %   p.d0_shallow = 0.1*p.d0_shallow;
end
% theta = 0*pi/180;
% p.ui=speed_steady*cos(theta); p.vi = speed_steady*sin(theta);

% Bottom profile
p.Dhole = 13e-3/p.xF;

if system == "channel"
    % geometry
    ch_width = 3;
    ch_length = 13;
    pillar_width = 1.5;
    pillar_pos = 8;
    ch_curve = 1;
    
    no_pillar = 0;
    %no_bottom = 0;
    
    %v0 = p.xF/p.TF
    
    channel_1 = (sqrt(p.yy.^2 + 0*p.xx.^2) < ch_width);
    channel_end = (p.xx > ch_length - 23);
    wall_cut = (max(abs(p.yy),abs(p.xx)) < 23.5) & (max(abs(p.yy),(p.xx)) < 18.5);
    channel_2 = (abs(p.yy ) - ch_width < ch_curve./abs(p.xx + 23 - ch_length));
    channel_2 = channel_2 | (no_bottom & (p.yy < 0 ));
    channel_3 = ( 3*abs(p.yy) + abs(p.xx + 23 - pillar_pos) ) < 4.5*ch_width;
    
    bdry = ( channel_1 | channel_2 | channel_end | channel_3 ) & wall_cut;
    
    
    delta = 0.3;
    point2=exp(-p.xx.^2/delta^2 - p.yy.^2/delta^2)/delta^2;
    reflective_bdry = conv2(bdry,point2,'same');
    normalize = max(max(reflective_bdry));
    reflective_bdry = fft2(1 - reflective_bdry/normalize);
    point2 = point2/normalize;
    deep = channel_1;
    
    pillar_1 = no_pillar | ( 3*abs(p.yy) + abs(p.xx + 23 - pillar_pos) ) > pillar_width;
    deep = ( ifft2(reflective_bdry)<1/2) & pillar_1;
    no_go_zone = (p.xx > 5) | (abs(p.yy) > 10);
elseif system == "slits"
    slit_width = 1.8;
    slit_x = 10;
    slit_spread = 2;
    slit_blockiness = 0.4;
    wall_cut = (max(abs(p.yy),abs(p.xx)) < 23.5) & (max(abs(p.yy),(p.xx)) < 20.5);

    slit_wall = abs(p.xx + 23 - slit_x)/slit_blockiness > 1;
    slit_wall = slit_wall | (abs(p.yy) - slit_spread/2 < slit_width);

    slit_middle = abs(p.xx + 23 - slit_x)/slit_blockiness > 1;
    slit_middle = slit_middle | (abs(p.yy) - slit_spread/2 > 0 );

    if no_bottom == 1
        slit_middle = slit_middle & ((abs(p.xx + 23 - slit_x)/slit_blockiness > 1) | (p.yy > 0));
    end

    bdry = slit_wall & wall_cut & slit_middle;

    delta = 0.3;
    point2=exp(-p.xx.^2/delta^2 - p.yy.^2/delta^2)/delta^2;

    reflective_bdry = conv2(bdry,point2,'same');
    normalize = max(max(reflective_bdry));
    reflective_bdry = fft2(1 - reflective_bdry/normalize);
    point2 = point2/normalize;
    deep = ( ifft2(reflective_bdry)<1/2);
    no_go_zone = (p.xx > 5) | (abs(p.yy) > 10);
end
%deep = 1 - ifft2(reflective_bdry);
%deep = ones(size(deep))
p.h = p.h0_deep.*deep + p.h0_shallow.*(1 - deep);
p.d = p.d0_deep.*deep + p.d0_shallow.*(1 - deep);
p.a = p.a0_deep.*deep + p.a0_shallow.*(1 - deep);



%% Execute

no_go_zone = fft2(conv2(no_go_zone,point2,'same'));

if p.useGPU == 1
    p.d = gpuArray(p.d);  p.a = gpuArray(p.a);
end

fname = [pwd,'/new-results-',convertStringsToChars(system),'/xi_',num2str(p.xi,'%.2f'),'_yi_',num2str(p.yi,'%.2f'),'_ui_',num2str(p.ui,'%.2f'),'_vi_',num2str(p.vi,'%.2f'),'_wall_',num2str(1 - no_bottom,'%.0f'),'_mem_',num2str(mem,'%.3f'),'.mat'];
file_exist =  exist(fname,'file');


if file_exist == 2 %File does not exists
    disp(['Result file already exists for input conditions. Returning.'])
    return
end

    [x_data,y_data,t_data, eta_data]  = trajectory(p,plotoption, width, type, no_go_zone); 
    %theta = atan(diff(x_data)./diff(y_data))*180/pi;  

    % Save data
    p = rmfield(p,'I_hat');
    p = rmfield(p,'K2');
    p = rmfield(p,'Kx');
    p = rmfield(p,'KxiKy');
    p = rmfield(p,'KxmiKy');
    p = rmfield(p,'Ky');
    p = rmfield(p,'abs_K');
    p = rmfield(p,'dissMatrix');
    p = rmfield(p,'dissMatrix_half');
    p = rmfield(p,'eta0');
  %  p = rmfield(p,'kx');
  %  p = rmfield(p,'ky');
    p = rmfield(p,'phi0');
    p = rmfield(p,'shift1');
    p = rmfield(p,'shift2');
    p = rmfield(p,'x');
    p = rmfield(p,'xxx');
    p = rmfield(p,'y');
    p = rmfield(p,'yyy');

    p = rmfield(p,'a');
    p = rmfield(p,'h');

    % save(fname,...
    %     'x_data','y_data','t_data','p','eta_data','theta');
    save(fname,...
        'x_data','y_data','t_data','p','eta_data');

   
end





