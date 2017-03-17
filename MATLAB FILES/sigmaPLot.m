clc;
% clear all;
close all;
cd ..;
fid = fopen('./OUTPUT/sigma.txt');
sigma = dlmread('./OUTPUT/sigma.txt','\t');
fclose(fid);


x = sigma(:,1)'./1e-3; %converting axis into meV
y = sigma(:,2)';
plot(x,y,'LineWidth',2);
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('Energy(meV)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('\sigma_{xy}(e^2/h)','FontSize',12,'FontWeight','bold','Color','k')
title('j/t = 3.5','FontSize',12,'FontWeight','bold','Color','k')
set(gca,'FontWeight','bold')