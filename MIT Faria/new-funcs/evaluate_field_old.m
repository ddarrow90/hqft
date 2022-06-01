function val = evaluate_field_old(eta,yi,xi,p)

if exist("p.Ky")
    Ky = p.Ky;
    Kx = p.Kx;
else
    kx  =  2*pi*1i/p.Lx*[0:p.Nx/2-1 0 -p.Nx/2+1:-1];
    ky  =  2*pi*1i/p.Ly*[0:p.Ny/2-1 0 -p.Ny/2+1:-1];
    Kx  = zeros(p.Ny,p.Nx); Ky = zeros(p.Ny,p.Nx);

    for i=1:p.Ny
        Kx(i,:) = kx;
    end
    
    for i=1:p.Nx
        Ky(:,i) = ky;
    end
end

ix  = find(p.xx(1,:)>xi,1); iy = find(p.yy(:,1)>yi,1);
if isempty(ix); ix=p.Nx; end
if isempty(iy); iy=p.Ny; end
shiftx = p.xx(1,ix)-xi;
shifty = p.yy(iy,1)-yi;

temp = real(ifft2(exp(-Kx.*shiftx-Ky.*shifty).*eta));
val = temp(iy,ix);
end