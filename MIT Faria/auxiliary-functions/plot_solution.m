function plot_solution(eta,xi,yi,t,p,cplot, history)
    figure(1)
    warning off;
    if p.useGPU ==1
        eta = gather(eta);
    end
    
    if nargin < 6
        cplot=[-0.003,0.003];
    end
    eta_interp = interp2(p.xx,p.yy,eta,p.xxx,p.yyy,'linear',0);
    h = pcolor(p.xxx,p.yyy,eta_interp);set(h,'edgecolor','none','FaceAlpha',0.95); grid off; 
    colorbar; caxis(cplot)
    title(['t=',num2str(t)]); hold on;
    if ~exist('history','var') || p.num_drops > 1
        plot(xi,yi,'k.','MarkerSize',10)
    else
        plot(history(:,1), history(:,2), 'k')
        plot(xi,yi,'k.','MarkerSize',10)
    end
    v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
    axis square;
    drawnow;  hold off; 
    colormap parula;
    
if p.plotPS == 1
    figure(2)
    clf
    PS = abs(fft2(eta));
    [C,h2] = contourf(imag(p.Kx)/(2*pi),imag(p.Ky)/(2*pi),PS/max(max(PS)));
    set(h2,'LineStyle','none');
    colorbar;
    axis square;
    drawnow;
end
    
    