clc;
cd ..;



%INPUT format for Berry Curvature calculation should be
%Kx Ky BerryCurvatureforBand1 BerryCurvatureforBand2 ......

%BerryCurvature  =dlmread('BerryCurvature.txt','\t');

fid = fopen('./OUTPUT/BerryCurvature.txt');
BerryCurvature = dlmread('./OUTPUT/BerryCurvature.txt','\t');


nx = 61; % number of site along kx
ny = 61; % number of site along ky
delta = pi/30; % intermediate distance between sites

x= 0:delta:(nx-1)*delta;
y= 0:delta:(ny-1)*delta;


   
e = 1;
z1 = zeros(nx,ny);
z2 = zeros(nx,ny);    


n=  4;  %n determines the nth band for which Berry Curvature would be calculated



 velocityField1(:,1) = BerryCurvature(:,n+2);



for i=1:nx:(nx*(ny-1)+1)
    
     h1 = velocityField1((i:(i+(nx-1))),1)';
     i22(1,e) = trapz(x,h1);
     z1(e,:) = h1;
     
         
    e = e+1;
end


chernNumber = trapz(y,i22)./(5*15e-10)^2/(2*pi)  %for 11 by 11 grid, (11*15e-10)^2
surf(x,y,z1);
xlabel('KX (A^{0})');
ylabel('KY (A^{0})');
s= strcat('Berry Curvature: Chern number for this band = ',num2str(chernNumber));
title(s);
colorbar;
grid on;
















