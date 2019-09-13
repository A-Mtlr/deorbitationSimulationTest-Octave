% Equation resolution of Newton's second law, ie second order nonlinear differential equation
% First approach with the Euler method

% But Euler's method is a first order method, hence the need to apply it twice to our equation

% Solve r1'(t) = r2(t) with r1(t0) = Rearth + r0 = 6371.000,0 + 400.000,0
% Solve r2'(t) = 1/r1(t).r2(t)² - GM/r1(t)² + C1.G(r1(t)) - C2.r2² with r2(t0) = 0
r1_t0 = 6771000;                  
r2_t0 = 0;      % radial velocity at 400km altitude

h = 0.1;                   % Time step
t = 0:h:500;               

distance = zeros(1,length(t));

%initialization of the list containing the distance to the center of the Earth at each instant t
r1star = zeros(1,length(t));  % Preallocate array
%initialization of the list containing the radial velocity at each instant t
r2star = zeros(1,length(t));
%initialization of the list containing the centrifugal force at each instant t
Fcentrifuge = zeros(1,length(t));
%initialization of the list containing the gravitational force at each moment t
Fgravity = zeros(1,length(t));
%initialization of the list containing the solar radiation force at each instant t
FsolarRadiation = zeros(1,length(t));
%initialization of the list containing the atmospheric drag force at each instant t
FdragAtm = zeros(1,length(t));


%comparative functions
expon = zeros(1,length(t));
logari = zeros(1,length(t));
lineaR = zeros(1,length(t));


r1star(1) = r1_t0;   % initial distance        
r2star(1) = r2_t0;   % initial velocity

% preliminary calculations 
GM=6.674*10^(-11)*5.972*10^24;
Pmf=0.0000034390458215;
d=0.01;
m=1.33;


% Solve r1 and r2 as respectively a first order linear and non linear differential equation
for i=1:(length(t))
    
    if (r1star(i)<=6371000+400000) && (r1star(i)>6371000)
        
        %First use Euler's method
        k1 = r2star(i);  % approx of the velocity at instant i-1
        r1star(i+1) = r1star(i) + k1*h; % Approximate solution for next value of y
        
        % preliminary calculations for solar radiation force
        % first parameter C1
        alpha=atan(d/(r1star(i)));
        theta=2*asin(6371000/(r1star(i)));
        thetaPrime = alpha + theta;
        A=2*d*(1/(tan(thetaPrime))+1)*cos(pi/2-thetaPrime);
        A1=(3+2*sqrt(2)+A/d)*d*d/2; 
        nu=1-theta/(2*pi);
        C1=Pmf*A1*nu/m;
        % second parameter eta
        thetaSecond=(5*pi-3*theta)/32;
        partieSup=cos(thetaSecond)-1.99*cos(asin((sin(thetaSecond))/1.99));
        partieInf=cos(thetaSecond)+1.99*cos(asin((sin(thetaSecond))/1.99));
        eta=partieSup/partieInf;


        % preliminary calculations for atm drag force
        %first parameter C2
        % C2=Cd*A2/m with Cd=0.80 the drag coefficient and A2=d*d*sqrt(2) the cross sectional area
        C2=1.131370849898476/m;
        % second parameter rho the density of the atm fluid, but we have not an expression of it for such altitude
        % thus, we use data base MSISE-90
        % The MSISE-90 model describes the densities in Earth's atmosphere from ground to thermospheric heights (~600km)
        % and it's made by NASA Goddard Space Flight Center
        % With this database, we draw a line between each data to get a density for each distance (even rational number)
        %
        %if our distance of Earth's center is between:
        %    x1=Rearth+400km 
        %and x2=Rearth+380km
        %with y1=atm density at x1 altitude=0.00000000000389 kg.m^-3
        %and  y2=atm density at x2 altitude=0.00000000000555 kg.m^-3
        if (r1star(i) <= 6371000+400000) && (r1star(i) > 6371000+380000)
          %slope=(y1-y2)/(x1-x2)
          slope=(0.00000000000389-0.00000000000555)/(400000-380000);
          %yIntercept=y1-slope*x1 (=y2-slope*x2)
          yIntercept=0.00000000000389-slope*(6371000+400000);
          %equation of a line to obtain the value of rho with only 2 other values
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%    
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+380km 
        %and x2=Rearth+360km
        %with y1=atm density at x1 altitude=0.00000000000555 kg.m^-3
        %and  y2=atm density at x2 altitude=0.00000000000799 kg.m^-3  
        elseif (r1star(i) <= 6371000+380000) && (r1star(i) > 6371000+360000)

          slope=(0.00000000000555-0.00000000000799)/(380000-360000);
          yIntercept=0.00000000000555-slope*(6371000+380000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%      
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+360km 
        %and x2=Rearth+340km
        %with y1=atm density at x1 altitude=0.00000000000799 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000000000116 kg.m^-3  
        elseif (r1star(i) <= 6371000+360000) && (r1star(i) > 6371000+340000)

          slope=(0.00000000000799-0.0000000000116)/(360000-340000);
          yIntercept=0.00000000000799-slope*(6371000+360000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%              
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+340km 
        %and x2=Rearth+320km
        %with y1=atm density at x1 altitude=0.0000000000116 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000000000172 kg.m^-3  
        elseif (r1star(i) <= 6371000+340000) && (r1star(i) > 6371000+320000)

          slope=(0.0000000000116-0.0000000000172)/(340000-320000);
          yIntercept=0.0000000000116-slope*(6371000+340000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                      
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+320km 
        %and x2=Rearth+300km
        %with y1=atm density at x1 altitude=0.0000000000172 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000000000258 kg.m^-3  
        elseif (r1star(i) <= 6371000+320000) && (r1star(i) > 6371000+300000)

          slope=(0.0000000000172-0.0000000000258)/(320000-300000);
          yIntercept=0.0000000000172-slope*(6371000+320000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                              
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+300km 
        %and x2=Rearth+280km
        %with y1=atm density at x1 altitude=0.0000000000258 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000000000394 kg.m^-3  
        elseif (r1star(i) <= 6371000+300000) && (r1star(i) > 6371000+280000)

          slope=(0.0000000000258-0.0000000000394)/(300000-280000);
          yIntercept=0.0000000000258-slope*(6371000+300000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                     
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+280km 
        %and x2=Rearth+260km
        %with y1=atm density at x1 altitude=0.0000000000394 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000000000616 kg.m^-3  
        elseif (r1star(i) <= 6371000+280000) && (r1star(i) > 6371000+260000)

          slope=(0.0000000000394-0.0000000000616)/(280000-260000);
          yIntercept=0.0000000000394-slope*(6371000+280000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                            
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+260km 
        %and x2=Rearth+240km
        %with y1=atm density at x1 altitude=0.0000000000616 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000000000991 kg.m^-3  
        elseif (r1star(i) <= 6371000+260000) && (r1star(i) > 6371000+240000)

          slope=(0.0000000000616-0.0000000000991)/(260000-240000);
          yIntercept=0.0000000000616-slope*(6371000+260000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+240km 
        %and x2=Rearth+220km
        %with y1=atm density at x1 altitude=0.0000000000991 kg.m^-3
        %and  y2=atm density at x2 altitude=0.000000000166 kg.m^-3  
        elseif (r1star(i) <= 6371000+240000) && (r1star(i) > 6371000+220000)

          slope=(0.0000000000991-0.000000000166)/(240000-220000);
          yIntercept=0.0000000000991-slope*(6371000+240000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+220km 
        %and x2=Rearth+200km
        %with y1=atm density at x1 altitude=0.000000000166 kg.m^-3
        %and  y2=atm density at x2 altitude=0.000000000291 kg.m^-3  
        elseif (r1star(i) <= 6371000+220000) && (r1star(i) > 6371000+200000)

          slope=(0.000000000166-0.000000000291)/(220000-200000);
          yIntercept=0.000000000166-slope*(6371000+220000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+200km 
        %and x2=Rearth+180km
        %with y1=atm density at x1 altitude=0.000000000291 kg.m^-3
        %and  y2=atm density at x2 altitude=0.000000000551 kg.m^-3  
        elseif (r1star(i) <= 6371000+200000) && (r1star(i) > 6371000+180000)

          slope=(0.000000000291-0.000000000551)/(200000-180000);
          yIntercept=0.000000000291-slope*(6371000+200000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+180km 
        %and x2=Rearth+160km
        %with y1=atm density at x1 altitude=0.000000000551 kg.m^-3
        %and  y2=atm density at x2 altitude=0.00000000118 kg.m^-3  
        elseif (r1star(i) <= 6371000+180000) && (r1star(i) > 6371000+160000)

          slope=(0.000000000551-0.00000000118)/(180000-160000);
          yIntercept=0.000000000551-slope*(6371000+180000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+160km 
        %and x2=Rearth+140km
        %with y1=atm density at x1 altitude=0.00000000118 kg.m^-3
        %and  y2=atm density at x2 altitude=0.00000000326 kg.m^-3  
        elseif (r1star(i) <= 6371000+160000) && (r1star(i) > 6371000+140000)

          slope=(0.00000000118-0.00000000326)/(160000-140000);
          yIntercept=0.00000000118-slope*(6371000+160000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+140km 
        %and x2=Rearth+120km
        %with y1=atm density at x1 altitude=0.00000000326 kg.m^-3
        %and  y2=atm density at x2 altitude=0.000000018 kg.m^-3  
        elseif (r1star(i) <= 6371000+140000) && (r1star(i) > 6371000+120000)

          slope=(0.00000000326-0.000000018)/(140000-120000);
          yIntercept=0.00000000326-slope*(6371000+140000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+120km 
        %and x2=Rearth+100km
        %with y1=atm density at x1 altitude=0.000000018 kg.m^-3
        %and  y2=atm density at x2 altitude=0.000000508 kg.m^-3  
        elseif (r1star(i) <= 6371000+120000) && (r1star(i) > 6371000+100000)

          slope=(0.000000018-0.000000508)/(120000-100000);
          yIntercept=0.000000018-slope*(6371000+120000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+100km 
        %and x2=Rearth+80km
        %with y1=atm density at x1 altitude=0.000000508 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0000168 kg.m^-3  
        elseif (r1star(i) <= 6371000+100000) && (r1star(i) > 6371000+80000)

          slope=(0.000000508-0.0000168)/(100000-80000);
          yIntercept=0.000000508-slope*(6371000+100000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+80km 
        %and x2=Rearth+60km
        %with y1=atm density at x1 altitude=0.0000168 kg.m^-3
        %and  y2=atm density at x2 altitude=0.000331 kg.m^-3  
        elseif (r1star(i) <= 6371000+80000) && (r1star(i) > 6371000+60000)

          slope=(0.0000168-0.000331)/(80000-60000);
          yIntercept=0.0000168-slope*(6371000+80000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+60km 
        %and x2=Rearth+40km
        %with y1=atm density at x1 altitude=0.000331 kg.m^-3
        %and  y2=atm density at x2 altitude=0.00407 kg.m^-3  
        elseif (r1star(i) <= 6371000+60000) && (r1star(i) > 6371000+40000)

          slope=(0.000331-0.00407)/(60000-40000);
          yIntercept=0.000331-slope*(6371000+60000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+40km 
        %and x2=Rearth+20km
        %with y1=atm density at x1 altitude=0.00407 kg.m^-3
        %and  y2=atm density at x2 altitude=0.0949 kg.m^-3  
        elseif (r1star(i) <= 6371000+40000) && (r1star(i) > 6371000+20000)

          slope=(0.00407-0.0949)/(40000-20000);
          yIntercept=0.00407-slope*(6371000+40000);
          rho=slope*r1star(i)+yIntercept;
          
    %______________________________________________________________________________%                                                   
        %elseif our distance of Earth's center is between:
        %    x1=Rearth+20km 
        %and x2=Rearth+0km
        %with y1=atm density at x1 altitude=0.0949 kg.m^-3
        %and  y2=atm density at x2 altitude=1.17 kg.m^-3  
        elseif (r1star(i) <= 6371000+20000) && (r1star(i) > 6371000+0)

          slope=(0.0949-1.17)/(20000-0);
          yIntercept=0.0949-slope*(6371000+20000);
          rho=slope*r1star(i)+yIntercept;
          
        end

          % calculations of Forces
          Fcentrifuge(i) = r2star(i)*r2star(i)/r1star(i);
          Fgravity(i) = -GM/(r1star(i)*r1star(i));
          FsolarRadiation(i) = -C1*(1+eta);
          FdragAtm(i) = (1/2*rho*C2*r2star(i)*r2star(i))*cos(pi/4);
          
          %Second use of Euler's method
          k2 = Fcentrifuge(i)+Fgravity(i)+FsolarRadiation(i)+FdragAtm(i);
          r2star(i+1) = r2star(i) + k2*h;
          
          %Resizing the radial distance, in terms of altitude in km
          distance(i) = (r1star(i)-6371000)/1000;
          
          %comparative functions
          %expon(i) = round(exp(i)); %pb -> "warning: opengl_renderer: data values greater than float capacity."
          expon(i) = 400-(1 + (i)/11.2 + ((i/11.2))^2/2); %Taylor's development
          logari(i) = -log(i)*GM;
          lineaR(i) = i;
    end
end

% To uncomment when displaying the force comparison (to comment otherwise)
%plot(distance, Fcentrifuge, distance, Fgravity, distance, FsolarRadiation, distance, FdragAtm);
%legend('Fcentrifuge', 'Fgravity', 'FsolarRadiation', 'FdragAtm');

% To uncomment when displaying the comparison of the radial distance and f(t) = R_{init} - exp(t/10)(to comment otherwise)
%plot(t, distance, t, expon); %, t, logari, t, lineaR);
%legend('distance radiale du CubeSat', 'f(t) = R_{init} - exp(t/10)');

plot(t, distance);
legend('radial distance of CubeSat', 'f(t) = ln(t)*GM');
