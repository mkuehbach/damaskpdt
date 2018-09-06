clear all
N = 100;
X=rand(N,1);
Y=rand(N,1); %sin(X*10*pi)+randn(size(X))/3; 
data=[X,Y];

plot(data(:,1),data(:,2),'.')
min(data)
max(data)



% apply routine
[bandwidth,density,X,Y]=kde2d(data);
% plot the data and the density estimate
surf(X,Y,density,'LineStyle','none'), view([0,90])
colormap parula, hold on, alpha(1.0)
colorbar


% show data points
set(gca, 'color', 'blue');
plot(data(:,1),data(:,2),'w.','MarkerSize',5)