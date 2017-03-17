clc;
close all;
clear all;

% cd ..;
fid2 = fopen('./Dispersion.txt');
Dispersion = dlmread('./Dispersion.txt','\t');
fclose(fid2);

Dispersion(:,1) = Dispersion(:,1)./15; %converting result for 1A0
for i=2:1:9000
    Energy = Dispersion(:,i)'./1e-3; %converting energy to meV
    plot(Dispersion(:,1)',Energy);
    hold on;
    
end

ylabel('Energy (meV)','FontSize',12,'FontWeight','bold','Color','k')
xlabel('J_{H}/t','FontSize',12,'FontWeight','bold','Color','k')
% ylim([0 5]);
%title('\gamma = 5x10^{-10}','FontSize',12,'FontWeight','bold','Color','k')
set(gca,'FontWeight','bold')