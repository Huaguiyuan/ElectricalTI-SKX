N = 21;
a = 1;


[x,y] = meshgrid((0:1:N-1)*a,(0:1:N-1)*a);

Sq = vec2mat(BoundState(:,1),N);
%Sq = max(BoundState') - Sq;
contourf(x,y,Sq);
colorbar;


xlabel('X');
ylabel('Y');
title('Psi Squared Plot: 21x21 SKX, first hole(\Gamma), j/t = 0.2');
clear all;