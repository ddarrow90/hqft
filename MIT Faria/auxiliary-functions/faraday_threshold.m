function [is_stable,eta_max] = faraday_threshold(p)

is_stable = 1;
% Remove drop
p.drop_force{1} = @(x,y) 0;      
t=0;

% Set initial condition 
phi     = zeros(size(p.xx));
pert    = 10^(-5);
eta     = pert*( 0*exp(-5*(sqrt((p.xx-1).^2+p.yy.^2)).^2) + besselj(0,(2*pi)*sqrt(p.xx.^2+p.yy.^2)) ); 
xi = 0; yi=0; ui=0.0; vi=0.0; 

% Fourier transform (2d) of initial condition
phi_hat = fft2(phi);               
eta_hat = fft2(eta);

for n=1:p.nimpacts
   
    eta = real(ifft2(eta_hat));
    warning('off');
    %plot_solution(eta,xi,yi,t,p)
    
    if mod(n-1,1)==0
        eta_max=max(max(abs(eta)));
        disp(['impact=',num2str(n),' eta_max=',num2str(eta_max)]);
    end
    
        if eta_max > 10*pert
            disp('Probably above threshold. Stopping simulation...');
            plot_solution(eta,xi,yi,t,p,[-10*pert,10*pert]);
            is_stable = 0;
            break
        elseif eta_max < pert/10
            disp('Probably below threshold. Stopping simulation...');
            plot_solution(eta,xi,yi,t,p,[-10*pert,10*pert]);
            is_stable = 1;
            break
            
        end
        
    % Evolve wave between impacts
    [phi_hat, eta_hat]           = evolve_wave(phi_hat, eta_hat, t, p);
    
    t = t+p.impact_interval;

end






