clc;
cd ..;
fid = fopen('./OUTPUT/BerryCurvature.txt');
BerryCurvature = dlmread('./OUTPUT/BerryCurvature.txt','\t');
%INPUT format for Berry Curvature calculation should be
%Kx Ky BerryCurvatureforBand1 BerryCurvatureforBand2 ......
%BerryCurvature  =dlmread('BerryCurvature.txt','\t');

nx = sqrt(size(BerryCurvature(:,1),1)); % number of site along kx
ny = sqrt(size(BerryCurvature(:,1),1)); % number of site along ky
kx = reshape(BerryCurvature(:,1),[nx,ny]);
ky = reshape(BerryCurvature(:,2),[nx,ny]);
Nx = 5;%number of site in real space in x or y direction
a = 15e-10; % intersite distance between two neighbouring points
x= kx(:,1)';
y= ky(1,:);


e = 1;
z1 = zeros(nx,ny);
z2 = zeros(nx,ny);    

n=  3;  %n determines the nth band for which Berry Curvature would be calculated
BerryCurvature2 = reshape(BerryCurvature(:,n+2),[nx,ny]);


for i=1:nx
    i22(1,e) = trapz(x,BerryCurvature2(:,n));
    e = e+1;
end


chernNumber = trapz(y,i22)./(Nx*a)^2/(2*pi)  %for 11 by 11 grid, (11*15e-10)^2
surf(x,y,BerryCurvature2);
xlabel('KX (A^{0})');
ylabel('KY (A^{0})');
s= strcat('Berry Curvature: Chern number for this band = ',num2str(chernNumber));
title(s);
colorbar;
grid on;
















