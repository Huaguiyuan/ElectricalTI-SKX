clc;
clear all;
close all;

s = load('matlab.mat');
a = s.d(:,1);
b = s.d(:,2);
plot(a',b');
plt = Plot()
plt.XLabel = 'Energy (meV)';
plt.YLabel = '\sigma_{xy} (e^2/h)';
plt. XLim = [0 5];
plt. YLim = [-2 1];

