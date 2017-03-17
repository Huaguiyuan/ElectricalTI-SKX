clc;

N = 6;      %number of skyrmion
nx = 11;
ny = 11*N;  %skyrmion nanoribbon is being generated in y direction 

delta = 15; %distance in angstorm between two neighbouring site
x = 0:delta:(nx-1)*delta;
y=  0:delta:(ny-1)*delta;

z = 0;
JH = -0.04;
t = -1;
Phi = -1;

fid=fopen('input.txt','w');
fprintf(fid,'%12s %12s %12s %11s %11s %11s %11s %11s %11s\n','x','y','z','sx','sy','sz','JH','t','Phi');



countSite = 1;

for y= 0:delta:(ny-1)*delta
    for x= 0:delta:(nx-1)*delta
        
        if (countSite > 121)
            countSite = 1;
        end
        
        fprintf(fid,'%12.7f %12.7f %12.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f\n',x,y,z,SKXtexture(countSite,4),SKXtexture(countSite,5),SKXtexture(countSite,6),JH,t,Phi);
        
        
        countSite = countSite + 1;
    end
    
end











 fprintf(fid,'\n');
 fprintf(fid,'\n');
 fprintf(fid,'\n');
 fprintf(fid,'\n'); 


fclose(fid);