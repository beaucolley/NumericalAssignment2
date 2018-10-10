clear all
close all

M = csvread("out_shock.csv");

mVals = M(:,1);
theta = M(:,2);
B_L = M(:,3);
B_U = M(:,4);

figure(1)
title("\theta - \beta - M");
xlabel("\theta");
ylabel("\beta");
hold on;

i = 1;
j = 1;

while i<length(M)
    while mVals(j) == mVals(j+1) && j<length(M)
        j = j+1;
    end
    plot(theta(i:j),B_L(i:j),'bo-',theta(i:j),B_U(i:j),'ro-');
%     plot(theta(i:j),B_U(i:j));
    disp("i = "+i);
    disp("j = "+j);
    i = 1 + j;
    j = j + 1;
end
    