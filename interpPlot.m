clear all
close all

M = csvread("out_interp_plot.csv");
N = csvread("in_interp.csv",1);

res = 0.1   
X = res*(1:1:length(M));

figure(1)
hold on
title("Lagrange Interpolation"); 
p1 = plot(X,M(:,1));
p2 = plot(N(:,1),N(:,2),'x');
xlabel("x");
ylabel("f(x)");

legend([p1,p2],["Interpolation","Data Points"]);

% figure(2)
% hold on;
% title("Cublic Spline Interpolation"); 
% xlabel("x");
% ylabel("f(x)");
% p2 = plot(X,M(:,2));
% plot(N(:,1),N(:,2));