% This is a simple script for plotting the Rosenblatt PDF for a particular value of D.
%% Compute PDF (this takes some time to run)
D = 0.3;
x = linspace(-2,5,30);

M = 50;
N = 5;
pdf = RosenblattPDF( x,D,M,N);

%% Plot
% Fits a spline to use inbetween function evaluations
% (requires curve fitting toolbox)
xx = linspace(-2,5,1000);
yy = spline(x,pdf,xx);
plot(x,pdf,'.',xx,yy,'MarkerSize',20);