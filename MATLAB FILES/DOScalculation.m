clc;
close all;

%this code calculates Density of States for a given bandstructure

Bandstructure = EigenValueDispersion(:,3:244);
nkx = 61;      %total number of points along kx 
nky = 61;      %total number of points along ky
NE = 242;      %total number of Energy point at each K



Emin = 0.02;      % min energy level to calculate DOS, unit eV
Emax = 0.12;      % max energy level to calculate DOS, unit eV 
gamma = 1e-3;     %broadening factor for calculating Density of States

E = linspace(Emin,Emax,1e5);
DOS = zeros(1,length(E));


for i =1:1:(nkx*nky)
    for j=1:1:NE
        
        
        Ek = Bandstructure(i,j);               % Ek contains an energy point from the bandstructure
        DOS = DOS +    (gamma/pi)./((E-Ek).^2+(gamma/2).^2); %Density of states being calculated using methods described in Applied QUANTUM mechanics, Levi. section 5.7.3
        
        
        
        
    end
end


plot(DOS,E);
xlabel('Density of State');
ylabel('Energy(eV)');