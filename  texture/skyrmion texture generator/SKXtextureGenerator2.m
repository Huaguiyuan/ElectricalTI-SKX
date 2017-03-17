clc;
clear all;
close all;

%Author: Tonmoy Kumar Bhowmick, Date: 12 February, 2016
%This code generates skyrmion texture from analytic expression of skyrmion
%Input has to be radius(lambda) of skyrmion, distance between grid points,
% total number of grid point in a square lattice,
%Nagaosa, Naoto, and Yoshinori Tokura. "Topological properties and dynamics
%of magnetic skyrmions." Nature nanotechnology 8.12 (2013): 899-911.


fid = fopen('input.txt','w'); %where Skyrmion texture would be stored

N = 11;  %total number of grid points in a square lattice, here N= 5 means 
        %5 by 5 square grid
        
        
a = 15; %intersite distance between grid points in angstorm, although unit
        %unit does not matter
        
        
lambda = ((N-1)/2)*a; %Skyrmion radius; assuming N would always be odd, skyrmion
                      %center would always be at zero/center
lambda2 = ((N-1)/2); %Skyrmion radius (in terms of grid points); assuming N would always be odd, skyrmion
                      %center would always be at zero/center
                      
m = 1;        %Different configuration of skyrmion
gamma = pi/2; %Different configuration of skyrmion
z = 0;        %vertical coordinate of skyrmion plane


for j= -((N-1)/2):1:((N-1)/2)
    for i= -((N-1)/2):1:((N-1)/2)
        
        phi = m*atan2(j,i) + gamma;
        r = sqrt(i.^2+j.^2);
        theta = pi*(1-r/lambda2);
        
        sx = cos(phi)*sin(theta);
        sy = sin(phi)*sin(theta);
        sz = cos(theta);
        fprintf(fid,'%8.3f  %8.3f  %8.3f %8.3f %8.3f %8.3f\n',(i+((N-1)/2))*a,(j+((N-1)/2))*a,z,sx,sy,sz);
        %shifting the skyrmion center coordinate to follow convention in
        %the input fule structure, the convention is the file will start
        %from (0,0) and x coordinate would continue to increase keeping y
        %fixed to generate the neighbour intex correctly in C++ code
    end
    
    
end

fclose(fid);

        


















