function [Fx,Fy] = compute_slope_eta(eta_hat,xi,yi,p)

% Surface gradient
    eta_x_hat = eta_hat*p.Kx_new;
    eta_y_hat = p.Ky_new*eta_hat;

% Shift fourier spectrum to nearest point
    for n=1:p.num_drops
        %{
        ix  = find(p.xx(1,:)>xi(n),1); iy = find(p.yy(:,1)>yi(n),1);
        if isempty(ix); ix=p.Nx; end
        if isempty(iy); iy=p.Ny; end
        shiftx = p.xx(1,ix)-xi(n);
        shifty = p.yy(iy,1)-yi(n);

        eta_x = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_x_hat)); % Shifts solution 
        eta_x1, temp = evaluate_field(eta_x_hat,yi(n),xi(n),p);
        Fx(n) = eta_x(iy,ix); 
        eta_y = real(ifft2(exp(-p.Kx.*shiftx-p.Ky.*shifty).*eta_y_hat));
        eta_y1 = evaluate_field(eta_y_hat,yi(n),xi(n),p);
        Fy(n) = eta_y(iy,ix);
        %}
        Fx(n) = evaluate_field(eta_x_hat,yi(n),xi(n),p);
        Fy(n) = evaluate_field(eta_y_hat,yi(n),xi(n),p);
    end
        
end