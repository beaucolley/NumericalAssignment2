clear all
close all

M = csvread("out_waveeqn_EU.csv");

dt = 10/length(M);
X = 1:1:length(M(1,:));
plotRes = 2;
steps = round(plotRes/dt);

figure(1);

for i=1:1:length(M)
   plot(X,M(i,:));
   pause(0.01);
end