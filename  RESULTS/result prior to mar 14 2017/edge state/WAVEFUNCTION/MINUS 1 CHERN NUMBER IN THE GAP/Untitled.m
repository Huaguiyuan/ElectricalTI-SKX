clc;
close all;
I=1:150;
plot(I,e(1,:),'r',I,flip(e(1,:)),'k');
% hold on;
% plot(flip(e(1,:)));
plt = Plot();
plt.XLabel = 'site number';
plt.YLabel = '|\psi|^2';