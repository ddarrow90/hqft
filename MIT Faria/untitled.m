rng(1)
xi = - 22 + 0.5*randn(10,1);
yi = 0 + 0.06*randn(10,1);
ui = 0.7 + 0.06*randn(10,1);
vi = 0.18 + 0.02*randn(10,1);


for i = 1:10
    hole(xi(i),yi(i),ui(i),vi(i),0,1);
    hole(xi(i),yi(i),ui(i),vi(i),0,0);
end


%convert_ICs_to_lshallow = 0;
for i = 1:10
    hole(xi(i),-yi(i),ui(i),-vi(i),0,0);
end

%% Plot 1
for i = 1:10
    fname = [pwd,'/new-results/xi_',num2str(xi(i),'%.2f'),'_yi_',num2str(yi(i),'%.2f'),'_ui_',num2str(ui(i),'%.2f'),'_vi_',num2str(vi(i),'%.2f'),'_wall_1.mat'];
    hold on;
    S(i) = load(fname);
    plot(S(i).x_data,S(i).y_data);
end

p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');

axis square;
    drawnow;  hold off; 

%% Plot 2
for i = 1:10
    fname = [pwd,'/new-results/xi_',num2str(xi(i),'%.2f'),'_yi_',num2str(yi(i),'%.2f'),'_ui_',num2str(ui(i),'%.2f'),'_vi_',num2str(vi(i),'%.2f'),'_wall_0.mat'];
    hold on;
    S(i) = load(fname);
    plot(S(i).x_data,S(i).y_data);
end

p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
    drawnow;  hold off;     