close all
clear all
clc
%%
i=[0.5751561E-01 28.24824 89.86656 97.39682 179.9425];
headANG=[0 30 93.5 101 180];
plot(i,headANG)
grid on
hold on
[p,S] = polyfit(i,headANG,3);
iMAX=[0 28 90 97.5 180];
[y_fit,delta] = polyval(p,i,S)
plot(i,y_fit,'r-')
plot(i,y_fit+2*delta,'m--',i,y_fit-2*delta,'m--')
title('Linear Fit of Data with 95% Prediction Interval')
legend('Data','Fit','95% Prediction Interval')
