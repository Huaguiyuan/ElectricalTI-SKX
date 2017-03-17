clc;
cd ..;

fid1 = fopen('./OUTPUT/BerryCurvature.txt');
BerryCurvature = dlmread('./OUTPUT/BerryCurvature.txt','\t');
fid2 = fopen('./OUTPUT/Dispersion.txt');
EigenValueDispersion = dlmread('./OUTPUT/Dispersion.txt','\t');

nx = 61; % number of site along kx
ny = 61; % number of site along ky
delta = pi/30; % intermediate distance between sites

x= 0:delta:(nx-1)*delta;
y= 0:delta:(ny-1)*delta;

x11 = BerryCurvature(:,3:52);
s1 ='Berry curvature for Band';
s3='BerryCurvature';
counter = 0;

for n=1:1:50
   
        e = 1;
        z1 = zeros(nx,ny);
        z2 = zeros(nx,ny);    
        %n= 1;



         velocityField1(:,1) = x11(:,n);



        for i=1:nx:(nx*(ny-1)+1)

             h1 = velocityField1((i:(i+(nx-1))),1)';
             i22(1,e) = trapz(x,h1);
             z1(e,:) = h1;


            e = e+1;
        end
        
        chernNumber(n,1) = trapz(y,i22)./(5*15e-10)^2/(2*pi);  %there might be ambiguity here, for pi (2 surface) or 2*pi (single surface)
        
        s2 = int2str(n-1);
        s= strcat(s1,s2,':chern number: ',num2str(chernNumber(n,1)));

%         fig = figure;
%         surf(x,y,z1);
%         xlabel('KX (A^{0})');
%         ylabel('KY (A^{0})');
%         title(s);
%         colorbar;
%         grid on;
%         s4 = int2str(counter);
%         counter = counter +1;
%         s5=strcat(s3,s4);
%         print(fig,s5,'-dtiff');

        



end

fclose(fid1);
fclose(fid2);

close all;





