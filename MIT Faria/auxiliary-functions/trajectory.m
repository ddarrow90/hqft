function [x_data,y_data,t_data,eta_data] = trajectory(p,plotgap, width, type, no_go_zone)

function plot_sol(ifield,p,x1_0,y1_0,x1,y1)
    clf;
    %surf(p.xx,p.yy,abs(phi_hat))
    %hold on;
    
    %axis([-24 24 -24 24 -1 1])
    surf(p.xx,p.yy,real(ifft2(ifield)))
    hold on;
    %surf(p.xx,p.yy,1 + real(ifft2(p.Kx.*ifield)))
    stem3(x1,y1,2,'filled')
    stem3(x1_0,y1_0,2,'--m')
    axis([-24 -20 -2 2 -0.3 0.3])
    colormap summer
    shading interp
    pause(0.2)
end

% Set initial condition 
t=0;
phi     = p.phi0;
eta     = p.eta0; 
psi     = p.psi0;
xi = p.xi; yi=p.yi; ui=p.ui; vi=p.vi; 
spin = p.spin;

% Create vector for position
x_data = zeros(p.nimpacts,p.num_drops); y_data = x_data; t_data=zeros(p.nimpacts,1);

% Create vector for wavefield
if p.store_wavefield == 0; 
    eta_data = zeros(p.Nx,p.Ny);
else
    eta_data = zeros(p.Nx,p.Ny,p.nimpacts); 
end

% % Move wave field to GPU
if p.useGPU == 1
    eta = gpuArray(eta);
    phi = gpuArray(phi);
    psi = gpuArray(psi);
    x_data = gpuArray(x_data);
    y_data = gpuArray(y_data);
    t_data = gpuArray(t_data);
end

% Fourier transform (2d) of initial condition
phi_hat = fft2(phi);               
eta_hat = fft2(eta);
psi_hat = fft2(psi);
history = [xi, yi];

for n=1:p.nimpacts
    
    if plotgap == 1 || mod(n-1,plotgap)==0
        eta = real(ifft2(eta_hat));   
        eta_max=max(max(abs(eta)));
        if eta_max > 1
            disp('Probably above threshold. Stopping simulation...');
            plot_solution(eta,xi,yi,t,p, [-1e-3, 1e-3], history);
            break;
        end   
    end
    
    
    % Store position (and possibly wavefield)
    x_data(n,:) = xi; y_data(n,:) = yi; t_data(n) = t;
    if ismember(n,p.store_wavefield) == 1 
        eta_data(:,:,n) = gather(eta); 
    else
        eta_data = eta_data*(n-1)/n + eta*1/n;
    end
    
    

     % Drop impact
    [ui, vi, phi_hat, psi_hat] = drop_impact(xi,yi, ui, vi, phi_hat, eta_hat, psi_hat, spin, p, width);
    
    if plotgap == 1 || mod(n-1,plotgap)==0
        [ui_0, vi_0] = compute_slope_eta(eta_hat,xi,yi,p);
        ui_0 = - p.G.*ui_0./p.cf_impact;
        vi_0 = - p.G.*vi_0./p.cf_impact;
        xi_0 = xi;
        yi_0 = yi;
        %     disp(['impact=',num2str(n)]);
    end
    
    % Evolve drops between impacts
    [xi, yi, ui, vi] = evolve_drops(xi, yi, ui, vi, p);
    % Add noise (if included)
    if p.sig_noise ~= 0
        ui = ui + p.sig_noise.*randn; vi = vi + p.sig_noise.*randn;
    end
    
    % Evolve wave between impacts
    [phi_hat, eta_hat, psi_hat]           = evolve_wave_BDF4(phi_hat, eta_hat, psi_hat, t, p, type);
    
        % Diplay results
        % Plot
    history = [history; xi, yi];
    if plotgap == 1 || mod(n-1,plotgap)==0
        max_val = 0.5*max(max(abs(ifft2(eta_hat))));
        plot_solution(eta,xi,yi,t,p, [-1e-3, 1e-3], history);
        %plot_sol(eta_hat,p,xi_0,yi_0,xi,yi);

        ui_av = xi - x_data(n,:); vi_av = yi - y_data(n,:);
        disp(['impact=',num2str(n), newline, '    pos = (',num2str(xi),', ',num2str(yi),...
        ')',newline, '    v_new = (',num2str(ui_av),', ',num2str(vi_av),...
        ')',newline, '    v_for = (',num2str(ui_0),', ',num2str(vi_0),...
        ')',newline, '    eta_max=',num2str(eta_max),...
        '  |eta_boundary|=',num2str(max(abs(eta(1,:))))]);
        %     disp(['impact=',num2str(n)]);
    end
    
    t = t+p.impact_interval;

    % Stop when drop goes out of domain
    ri  = sqrt(xi.^2+yi.^2);
    ri0 = sqrt(x_data(1,:).^2+y_data(1,:).^2);
    
    if evaluate_field(no_go_zone, yi, xi, p) > 1/2
        disp('Stopping simulation.');
        plot_solution(eta,xi,yi,t,p, [-1e-3, 1e-3], history);
        x_data(n:end,:) = [];
        y_data(n:end,:) = [];
        t_data(n:end,:) = [];   
        break
    end


end

if p.useGPU ==1
    x_data = gather(x_data); y_data = gather(y_data); 
    eta_data = gather(eta_data);
    t_data = gather(t_data);
end

end





