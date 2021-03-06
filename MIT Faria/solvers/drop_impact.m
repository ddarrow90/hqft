function [ui, vi, phi_hat, psi_hat] = drop_impact(xi,yi, ui, vi, phi_hat, eta_hat, psi_hat, spin, p, delta)
% DESCRIBES THE INSTANTANEOUS INTERACTION BETWEEN THE DROP AND THE SURFACE

[Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p);

% Drop speed immediately after impact
 ui = - p.G.*Fx./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*ui;
 vi = - p.G.*Fy./p.cf_impact.*(1-exp(-p.cf_impact)) + exp(-p.cf_impact).*vi;
% Wave immediately after impact
for n=1:p.num_drops
    %phi_hat =  phi_hat + p.M(n)*p.G/(p.hx*p.hy).*exp(-p.Kx*(p.Lx/2+xi(n))-p.Ky*(p.Ly/2+yi(n))); 
    point2=fft2(exp(-p.xx.^2/delta - p.yy.^2/delta));
    point2=exp(-(p.Kx)*(xi(n))-(p.Ky)*(yi(n))).*point2;
    phi_hat = phi_hat + p.M(n)*p.G/(p.hx*p.hy).*point2;
    % this is _fine_, but we can do better
end

end