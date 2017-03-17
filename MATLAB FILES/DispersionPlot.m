clc;
close all;
clear all;
cd ..;
fid = fopen('./OUTPUT/Dispersion.txt');
Dispersion = dlmread('./OUTPUT/Dispersion.txt','\t');



nx = sqrt(size(Dispersion(:,1),1)); % number of site along kx
ny = sqrt(size(Dispersion(:,1),1)); % number of site along ky
kx = reshape(Dispersion(:,1),[nx,ny]);
ky = reshape(Dispersion(:,2),[nx,ny]);

start_band = 24;  % band indexing start from 0 to 241 for total of 242 bands
end_band = 30;

for j= (start_band+3):1:(end_band+3)
        n= j;
        z1 = reshape(Dispersion(:,n),[nx,ny]);
        z1 = z1./1e-3;  %converting into meV
        surf(kx,ky,z1);

        hold on;

end
title('J/t = 0.05');
xlabel('Kx (A^{0})');
ylabel('ky (A^{0})');
zlabel('Energy(meV)');