function [X,Y] = fixed_potential(x,y,u,v, eta,T, dt,p)
% Evolves wave field in F-space through nsteps of size dt
X = zeros(max(size(x)),T);
Y = zeros(max(size(x)),T);
X(:,1) = x;
Y(:,1) = y;

for t=1:T-1
    for n=1:max(size(x))
        [Fx,Fy] = compute_slope_eta(fft2(eta),X(n,t),Y(n,t),p);
        
        u(n) = u(n) - dt*Fx;
        X(n,t+1) = X(n,t) + dt*u(n);

        v(n) = v(n) - dt*Fy;
        Y(n,t+1) = Y(n,t) + dt*v(n);
    end
end      

    
end