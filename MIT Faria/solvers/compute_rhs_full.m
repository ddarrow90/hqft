function [rhs1, rhs2, rhs3] = compute_rhs_full(phi_hat, eta_hat, psi_hat, t, p)
% Computes the RHS of the hyperbolic part
    
     phi_x = real(ifft2(p.Kx.*phi_hat));
     phi_y = real(ifft2(p.Ky.*phi_hat));
     psi_x = real(ifft2(p.Kx.*psi_hat));
     psi_y = real(ifft2(p.Ky.*psi_hat));
     
     u = phi_x + psi_y;
     v = phi_y - psi_x;

scaling = 10;
function out = mat_drv(eta_hat)
    out = zeros(p.Ny,p.Nx);
    %out = ifft2(p.Kx.*eta_hat).*u + ifft2(p.Ky.*eta_hat).*v;
    %out = fft2(out)/scaling;
end



rhs1 = - mat_drv(phi_hat) - p.g(t).*eta_hat + p.Bo*p.K2.*eta_hat + 2*p.nu0.*p.K2.*phi_hat;

rhs2 = - mat_drv(eta_hat) + DtN(phi_hat,p) + 2*p.nu0.*p.K2.*eta_hat; % DtN computes phi_z_hat  

rhs3 = zeros(p.Ny,p.Nx);%- mat_drv(psi_hat) + p.nu0.*p.K2.*psi_hat;
end