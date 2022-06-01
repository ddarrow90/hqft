function [rhs2] =  DtN_new(phi_hat,p)
% Approximation to \phi_z (i.e. Dirichelt-to-Neumann operator)
    %rhs2 = -p.Kx.*fft2(p.d.*ifft2(p.Kx.*phi_hat));
    %rhs2 = rhs2 - p.Ky.*fft2(p.d.*ifft2(p.Ky.*phi_hat));
    rhs2 = fft2(p.d.*ifft2(p.abs_K.*phi_hat));

end