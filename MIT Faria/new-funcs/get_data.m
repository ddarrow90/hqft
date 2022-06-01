rng(1)
xi = - 21 + 0.5*randn(10,1);
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


%% Plot 2
for i = 1:10
    fname = [pwd,'/new-results - Copy/xi_',num2str(xi(i),'%.2f'),'_yi_',num2str(yi(i),'%.2f'),'_ui_',num2str(ui(i),'%.2f'),'_vi_',num2str(vi(i),'%.2f'),'_wall_0.mat'];
    hold on;
    S(i) = load(fname);
    color{i} = [rand(),rand(),rand()];
    plot(S(i).x_data,S(i).y_data, 'Color', color{i});
    
end

p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
    drawnow;  hold off;     

figure;

for i = 1:10
    fname = [pwd,'/new-results - Copy/xi_',num2str(xi(i),'%.2f'),'_yi_',num2str(yi(i),'%.2f'),'_ui_',num2str(ui(i),'%.2f'),'_vi_',num2str(vi(i),'%.2f'),'_wall_1.mat'];
    hold on;
    S(i) = load(fname);
    plot(S(i).x_data,S(i).y_data,'Color',color{i});
end

p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
    drawnow;  hold off;     

%% Test 2

y1i = -1.5:0.2:1.5;
x1i = - 21.5*ones(size(y1i));
u1i = 0.7*ones(size(y1i));
v1i = 0*ones(size(y1i));

for i = 1:length(y1i)
    hole(x1i(i),y1i(i),u1i(i),v1i(i),0,0.905,"channel");
end

%% Plot 3
for i = 1:length(y1i)
    fname = [pwd,'/new-results-channel/xi_',num2str(x1i(i),'%.2f'),'_yi_',num2str(y1i(i),'%.2f'),'_ui_',num2str(u1i(i),'%.2f'),'_vi_',num2str(v1i(i),'%.2f'),'_wall_1.mat'];
    hold on;
    S(i) = load(fname);
    plot(S(i).x_data,S(i).y_data, 'Color',[i/length(y1i),0.3, 1 - i/length(y1i)], 'LineWidth', 1);
end

p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
    drawnow;  hold off;    

%% Plot 4
i0 = length(y1i)/2 + 1/2;
 sig = 100;
eta_new = zeros(p.Ny,p.Nx);
for i = 1:length(y1i)
    fname = [pwd,'/new-results-channel/xi_',num2str(x1i(i),'%.2f'),'_yi_',num2str(y1i(i),'%.2f'),'_ui_',num2str(u1i(i),'%.2f'),'_vi_',num2str(v1i(i),'%.2f'),'_wall_1.mat'];
    hold on;
    S(i) = load(fname);
    eta_new = eta_new + exp(-(i-i0)^2/sig)*(S(i).eta_data);
end
h = pcolor(p.xx,p.yy,eta_new);set(h,'edgecolor','none','FaceAlpha',0.95); grid off; 
colorbar; caxis([0 norm(abs(eta_new))/100])
%surf(p.xx,p.yy,eta_new);
p = S(1).p;
%v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    %contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
%colormap summer
%shading interp
axis square;
%    drawnow;  hold off;    
for i = 1:length(y1i)
    fname = [pwd,'/new-results-channel/xi_',num2str(x1i(i),'%.2f'),'_yi_',num2str(y1i(i),'%.2f'),'_ui_',num2str(u1i(i),'%.2f'),'_vi_',num2str(v1i(i),'%.2f'),'_wall_1.mat'];
    hold on;
    S(i) = load(fname);
    plot(S(i).x_data,S(i).y_data, 'Color',[1,0.3, 0], 'LineWidth', 1);
end

p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
    drawnow;  hold off;    
%% Fixed Potential
delta = 0.01;
point2=exp(-p.xx.^2/delta^2 - p.yy.^2/delta^2)/delta^2;
potfield = conv2(eta_new,point2,'same');
potfield = potfield/norm(potfield);


%potfield =eta_new/norm(eta_new);
[X,Y] = fixed_potential(x1i,y1i,u1i,v1i,3.5*potfield,400,0.1,p);

for i=1:size(X,1)
    plot(X(i,:),Y(i,:),'Color',[i/length(y1i),0.3, 1 - i/length(y1i)], 'LineWidth', 1)
    hold on;
end
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
    drawnow;  hold off;   
%% Test 3

y1i = 1.5:-0.6:0.3;
x1i = - 21.5*ones(size(y1i));
u1i = 0.7*ones(size(y1i));
v1i = 0*ones(size(y1i));
mem = 0.88:0.02:0.92;
for m = flip(mem)
    for i = 1:length(y1i)
        hole(x1i(i),y1i(i),u1i(i),v1i(i),1,m, "channel");
    end
end

for m = mem
    for i = 1:length(y1i)
        hole(x1i(i),y1i(i),u1i(i),v1i(i),0,m, "channel");
    end
end


%% Plot 5

for m = mem
    figure;
    for i = 1:length(y1i)
        fname = [pwd,'/new-results-channel/xi_',num2str(x1i(i),'%.2f'),'_yi_',num2str(y1i(i),'%.2f'),'_ui_',num2str(u1i(i),'%.2f'),'_vi_',num2str(v1i(i),'%.2f'),'_wall_1','_mem_',num2str(m,'%.3f'),'.mat'];
        hold on;
        S(i) = load(fname);
        plot(S(i).x_data,S(i).y_data, 'Color',[i/length(y1i),0.3, 1 - i/length(y1i)], 'LineWidth', 1);
        title(['mem = ',num2str(m,3)])
    end
    
    p = S(1).p;
    v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
        contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
    axis square;
    drawnow;    
end


%% Test slits

y1i = 1:0.075:2.8;
x1i = - 21.5*ones(size(y1i));
u1i = 0.7*ones(size(y1i));
v1i = 0*ones(size(y1i));

for i = 1:length(y1i)
    hole(x1i(i),y1i(i),u1i(i),v1i(i),0, 0.905,"slits");
end

for i = 1:length(y1i)
%    hole(x1i(i),y1i(i),u1i(i),v1i(i),1, 0.905,"slits");
end

%% Plot slits
data = [];

 figure;
 i0 = 0;
 sig = 100;
if exist("eta_new")
    eta_new = 0*eta_new;
end
for i = 1:length(y1i)
    fname = [pwd,'/new-results-slits/xi_',num2str(x1i(i),'%.2f'),'_yi_',num2str(y1i(i),'%.2f'),'_ui_',num2str(u1i(i),'%.2f'),'_vi_',num2str(v1i(i),'%.2f'),'_wall_1','_mem_',num2str(0.905,'%.3f'),'.mat'];
    hold on;
    S(i) = load(fname);
    if ~exist("p")
        p = S(1).p;
        eta_new = zeros(p.Ny,p.Nx);
    end
    %plot(S(i).x_data,S(i).y_data, 'Color',[0.5 + 0.5*i/length(y1i),0.3, 0.5 - 0.5*i/length(y1i)], 'LineWidth', 1);
    %plot(S(i).x_data,-S(i).y_data, 'Color',[0.5 - 0.5*i/length(y1i),0.3, 0.5 + 0.5*i/length(y1i)], 'LineWidth', 1);
    data(2*i-1) = (S(i).y_data(end) - S(i).y_data(end-1))/(S(i).x_data(end) - S(i).x_data(end-1));
    data(2*i) = -(S(i).y_data(end) - S(i).y_data(end-1))/(S(i).x_data(end) - S(i).x_data(end-1));
    eta_new = eta_new + exp(-(i - i0)^2/sig)*S(i).eta_data;
    eta_new = eta_new + exp(-(i- i0)^2/sig)*flip(S(i).eta_data);
end
h = pcolor(p.xx,p.yy,eta_new);set(h,'edgecolor','none','FaceAlpha',0.95); grid off; 
colorbar; caxis([0 norm(abs(eta_new))/100])
for i = 1:length(y1i)
   % plot(S(i).x_data,S(i).y_data, 'r', 'LineWidth', 1);
   % plot(S(i).x_data,-S(i).y_data, 'r', 'LineWidth', 1);
end
p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
drawnow;  
figure;
hist((data), 30)
%% Plot shitty slits

data = [];
 figure;
 i0 = 2;
 sig = 500;
if exist("eta_new")
    eta_new = 0*eta_new;
end
for i = 1:length(y1i)
    fname = [pwd,'/new-results-slits/xi_',num2str(x1i(i),'%.2f'),'_yi_',num2str(y1i(i),'%.2f'),'_ui_',num2str(u1i(i),'%.2f'),'_vi_',num2str(v1i(i),'%.2f'),'_wall_0','_mem_',num2str(0.905,'%.3f'),'.mat'];
    hold on;
    S(i) = load(fname);
    if ~exist("p")
        p = S(1).p;
        eta_new = zeros(p.Ny,p.Nx);
    end
   % plot(S(i).x_data,S(i).y_data, 'Color',[0.5 + 0.5*i/length(y1i),0.3, 0.5 - 0.5*i/length(y1i)], 'LineWidth', 1);
   % plot(S(i).x_data,-S(i).y_data, 'Color',[0.5 - 0.5*i/length(y1i),0.3, 0.5 + 0.5*i/length(y1i)], 'LineWidth', 1);
    eta_new = eta_new + exp(-(i - i0)^2/sig)*S(i).eta_data;
    data(i) = (S(i).y_data(end) - S(i).y_data(end-1))/(S(i).x_data(end) - S(i).x_data(end-1));
end
%h = pcolor(p.xx,p.yy,eta_new);set(h,'edgecolor','none','FaceAlpha',0.95); grid off; 
%colorbar; caxis([0 norm(abs(eta_new))/100])
for i = 1:length(y1i)
    plot(S(i).x_data,S(i).y_data, 'r', 'LineWidth', 1);
   % plot(S(i).x_data,-S(i).y_data, 'r', 'LineWidth', 1);
end
p = S(1).p;
v=[(p.d0_shallow+p.d0_deep)*0.49,(p.d0_shallow+p.d0_deep)*0.51];
    contour(p.xx,p.yy,p.d,v,'LineWidth',2,'LineColor','k');
axis square;
drawnow;  
figure;
hist(atan(data), 10)

