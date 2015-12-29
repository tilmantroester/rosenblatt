D = 0.3
x = linspace(-2,5,30);

M = 50;
N = 5;
pdf = RosenblattCDF( x,D,M,N,'pdf');

%% Plot
% Fits a spline to use inbetween function evaluations
% (requires curve fitting toolbox)
xx = linspace(-2,5,1000);
yy = spline(x,pdf,xx);
plot(x,pdf,'.',xx,yy,'MarkerSize',20);