function val = evaluate_field(eta,yi,xi,p)
imaginary = 0;


%val = (exp(-p.ky*(p.Ly/2+yi))/p.Ny)*eta*(exp(-p.kx*(p.Lx/2+xi)).'/p.Nx);
if exist("p.ky")
    ky = p.ky;
    kx = p.kx;
else
    kx  =  2*pi*1i/p.Lx*[0:p.Nx/2-1 0 -p.Nx/2+1:-1];
    ky  =  2*pi*1i/p.Ly*[0:p.Ny/2-1 0 -p.Ny/2+1:-1];
end

val = (exp(-ky*(p.Ly/2-yi))/p.Ny)*eta*(exp(-kx*(p.Lx/2-xi)).'/p.Nx);
if (imaginary == 0)
    val = real(val);
end
%val = (exp(-p.Ky.'*(p.Ly/2+yi))/p.Ny)*eta*(exp(-p.Kx.'*(p.Lx/2+xi))/p.Nx);
end