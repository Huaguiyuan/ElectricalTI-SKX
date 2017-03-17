%this code calculates sigma_xy (conductivity) for a given structure, lets
%say 11 1 lattice. the inputs are energy dispersion, chern number, berry
%curvature for a particular structure. Energy dispersion, Berry curvature
%for each band, chern number has to imported to the workspace to run this
%script. the variable velocityfield is meaningless in this respect.

clc;
cd ..;
fid = fopen( './OUTPUT/sigma.txt','wt');

nx = 5;  % total number of point along x direction
ny = 5;  % total number of point along y direction
basis = 2;  % for spin basis = 2
totalNumberOfBands = nx*ny*basis;



%%%%% start generaing k grid
nkx = 61; % number of site along kx of Brillouin zone
nky = 61; % number of site along ky of Brillouin zone
delta = pi/30; % intermediate distance between k point
kx= 0:delta:(nkx-1)*delta;
ky= 0:delta:(nky-1)*delta;
%%%%% end generating k grid

for i5= 25:26 %band indexes over which bands the sigma_xy will be calculated, band indexes are from 1 to 242

        nthBand = i5;  %bands are numbered from 1 2 ... to totalNumberOfBands (242)
        Emin = min(EigenValueDispersion(:,nthBand+2)); % EigenValueDispersion: kx ky band1 band2 band3 ...   bandN
        Emax = max(EigenValueDispersion(:,nthBand+2));
        E = linspace(Emin,Emax,10);
        offset = sum(chernNumber(1:(nthBand-1),1));
        
        
        for i4= E  %i4 contains the current value of energy for a paricular band above which the corresponding berry curvature for the same band should be set to zero

            EnergyDispersionOfnthBand = EigenValueDispersion(:,nthBand+2);
            EnergyLessthanparticularValue = (EnergyDispersionOfnthBand < i4); %determining which kx ky has to be selected such that corresponding energies fall below the energy range specified
            bc1 = BerryCurvature(:,nthBand+2);
            bc2 = bc1.*EnergyLessthanparticularValue; %bc2 contains the berry curvatures for which the integration needs to be done



            e = 1;
            z1 = zeros(nkx,nky);
            for i=1:nkx:(nkx*(nky-1)+1)

                 h1 = bc2((i:(i+(nkx-1))),1)';
                 i22(1,e) = trapz(kx,h1);
                 e = e+1;
            end

           chernNumber2 = trapz(ky,i22)./(5*15e-10)^2/(2*pi)+ offset; %for 11 by 11 grid, (11*15e-10)^2,  %there might be ambiguity here, for pi (2 surface) or 2*pi (single surface)
           fprintf( fid, '%f\t%f\n', i4,chernNumber2); 

        end
        
        bandgap = min(EigenValueDispersion(:,nthBand+3))- max(EigenValueDispersion(:,nthBand+2)) %determining band gap between nth and nth+1 band
        
        if (bandgap>(0.001))   %minimum bandgap has to be 1mev
            
            Emin2 = max(EigenValueDispersion(:,nthBand+2));
            Emax2 = min(EigenValueDispersion(:,nthBand+3));
            E2 = linspace(Emin2,Emax2,10);
            chernNumber3 = sum(chernNumber(1:(nthBand),1));
            
            for i6 = E2
                
                fprintf( fid, '%f \t%f\n', i6,chernNumber3);
                
 
            end
            
            
            
            
            
            
        end
        
        
        

end
fclose(fid);

