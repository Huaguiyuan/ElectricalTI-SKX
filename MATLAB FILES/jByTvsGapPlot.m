clc;
close all;
clear all;

cd ..;
fid2 = fopen('./OUTPUT/Dispersion.txt');
Dispersion = dlmread('./OUTPUT/Dispersion.txt','\t');
fclose(fid2);

Dispersion(:,1) = Dispersion(:,1).*1.5; %converting result for 1A0
for i=2:1:2970
    Energy = Dispersion(:,i)'./1e-3; %converting energy to meV
    plot(Dispersion(:,1)',Energy);
    hold on;
    
end

ylabel('Energy (meV)','FontSize',12,'FontWeight','bold','Color','k')
xlabel('J_{prime}','FontSize',12,'FontWeight','bold','Color','k')
ylim([0 15]);
title('SKX 5x5 supercell','FontSize',12,'FontWeight','bold','Color','k')
set(gca,'FontWeight','bold')