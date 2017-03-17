N =  5; %total number of sites in each direction

a= 1 ; %intersite distane in angstrom



[x,y] = meshgrid((0:1:N-1)*a,(0:1:N-1)*a);
SX = vec2mat(SpinTexture(:,1),N);
SY = vec2mat(SpinTexture(:,2),N);
SZ = vec2mat(SpinTexture(:,3),N);


contourf(x,y,SZ);
colorbar;
hold on;
quiver(x,y,SX,SY,'LineWidth',1.5);
xlim([0 N]);
ylim([0 N]);
xlabel('X');
ylabel('Y');
title('Spin Texture Plot: 5x5 TI(FM), j/t = 0.5, 1st e');
clear all;