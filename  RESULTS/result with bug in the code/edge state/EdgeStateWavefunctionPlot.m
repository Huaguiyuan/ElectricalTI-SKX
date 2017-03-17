clc;
close all;
    



for i=790:1:890
    current = x(:,i);
    b= reshape(current,[5,150]);
    figure;
    contourf(b);
    colorbar;
    
end