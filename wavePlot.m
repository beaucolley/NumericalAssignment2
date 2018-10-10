clear all
close all

M = csvread("out_waveeqn_EU.csv");
N = csvread("out_waveeqn_LW.csv");

dt = 10/length(M);
X = 1:1:length(M(1,:));
plotRes = 0.1;
plotSamples = 5
steps = round(plotRes/dt);

figure(1);
title("Euler Upwinding");
xlabel("x");
ylabel("f(x)");

figure(2);
title("Lax-Wendroff Method");
xlabel("x");
ylabel("f(x)");

for i=1:steps:5*steps
   figure(1);
   hold on;
   plot(X,M(i,:));
   figure(2);
   hold on
   plot(X,N(i,:));
   pause(0.01);
end