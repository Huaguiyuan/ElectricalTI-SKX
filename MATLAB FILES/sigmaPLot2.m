clc;
% clear all;
close all;
cd ..;
fid = fopen('./OUTPUT/sigma.txt');
sigma = dlmread('./OUTPUT/sigma.txt','\t');
fclose(fid);


[b117, I117] = sort(sigma(:,1));
c117 = sigma(:,2);
d117 = c117(I117);
plot(b117',d117');

x = b117'./1e-3;
y = d117';
plot(x,y,'LineWidth',2);
xlim([min(x) max(x)]);
ylim([min(y) max(y)]);
xlabel('Energy(meV)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('\sigma_{xy}(e^2/h)','FontSize',12,'FontWeight','bold','Color','k')
title('j/t = 0.05','FontSize',12,'FontWeight','bold','Color','k')
set(gca,'FontWeight','bold')