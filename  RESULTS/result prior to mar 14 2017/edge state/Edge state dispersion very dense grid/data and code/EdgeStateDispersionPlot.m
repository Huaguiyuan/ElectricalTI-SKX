
clc;
close all;
clear all;


% nx = 61; % number of site along kx
% 
% delta = pi/30; % intermediate distance between sites
% 
% x= 0:delta:(nx-1)*delta;

fid = fopen('EigenValueDispersion.txt');

EigenValueDispersion = dlmread('EigenValueDispersion.txt','\t');
[m1 n1] =  size(EigenValueDispersion);

start_band = 0;  % band indexing start from 0 to 241 for total of 242 bands
end_band = n1-4;

% j= start_band:1:end_band;

for j= (start_band+3):1:(end_band+3)
    
    
        n= j;
        EigenValueDispersion(:,n) = EigenValueDispersion(:,n)./1e-3; %converting energy(eV) to meV
        plot(EigenValueDispersion(:,1)',EigenValueDispersion(:,n)','r');
        hold on;

end

%title('Edge state for SKX on TI: velocity = 10%, J_{H}= -0.01, 5 skyrmion in nanoribbon');
xlim([0 6.2]);
ylim([0 6.2]);
ylabel('Energy(meV)','FontSize',18,'FontWeight','bold','Color','k')
xlabel('kx*a','FontSize',18,'FontWeight','bold','Color','k')
set(gca,'FontWeight','bold','fontsize', 15)


fclose(fid);










