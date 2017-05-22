%% Notes
%units based on data given -- mol/s, gram of catalyst, kPa
%check all variable values are reasonable; magnitude and sign
%% Component definitions
% 1 = cis-butene
% 2 = trans-butene
% 3 = 1-butene
% 4 = water  -- or w
% 5 = SBA    -- or a
% 6 = 1,3-butadiene
% 7 = i-butane
%% Inlet
F0 = [6.0230139,14.85743194,16.9150625,151.154161,0.00677639,1.827233333,0.18854306];  %mol/s
F0(4) = 10*sum(F0(1:3)); %Steam:Butene ratio
FT0 = sum(F0); %total molar flowrate in a given step
x0 = F0./FT0; %vector of compositions
%% Variable parameters
T0 = 491; %K
P0 = 7.5*10^2; %kPa
Height = 4; %m
TubeInnerDiameter = 0.50; %m
TubeThickness = 0.020; %m
NumberOfTubes = 5; %number of tubes
dz = 0.0001;   %m,height step size
Voidage = 0.4;    %voidage of packed bed
%% Reaction parameters
FlowArea = NumberOfTubes*pi*0.25*TubeInnerDiameter^2; %m2
k = 1.2*10^-7; %k is constant bc E=0 for all rxns
Kw0 = [3*10^-31, 1*10^-30, 1.2*10^-31]; %kPa^-3
dH = [-2.6*10^5,-2.2*10^5,-2.47*10^5]; %J/mol
KRatio = [150,80,100];
R = 8.3145; %J/mol K
CatalystRho = 2500; %kg/m3
ParticleDiameter = 0.001; %m
CatalystMass = CatalystRho*FlowArea*Height*(1-Voidage); %kg
N1 = Height/dz;   %number of steps
%% Reactor Model
Table = zeros;
Table(1,2) = T0;
Table(1,3) = P0;
Table(1,4) = x0(1);
Table(1,5) = x0(2);
Table(1,6) = x0(3);
Table(1,7) = x0(4);
Table(1,8) = x0(5);
for i = 1:N1
    %first iteration want inlet conditions
    if i == 1        
        T = T0; %K
        P = P0; %kPa
        F = F0; %mol/s
        CoolantTemp = 400; %K
    end
    %% Reaction
    FT = sum(F); %total molar flowrate in a given step
    x = F./FT; %vector of compositions
    pc4 = x(1:3).*P; %kPa
    pw = x(4)*P;     %kPa
    pa = x(5)*P;     %kPa
    Kw = Kw0.*exp(-dH./(R*T)); %kPa^-3

    rs = (k.*pc4.*pw^2)./((1+Kw.*pw^3).*(pw+KRatio.*pa)); %rate//per GRAM of catalyst
    m = 1000*CatalystMass/N1; %GRAMS of catalyst in each step
    r = rs.*m;   %number of moles of given butene consumed per second in a given step
 
    if any(r<0)
        errordlg('Negative rate; look at pressure drop')
        break
    end
    
    dF = [-r(1),-r(2),-r(3),-r(1)-r(2)-r(3),r(1)+r(2)+r(3),0,0]; %Assuming no side reactions
    %% Temperature Change
    DH = r.*(-5*10^4); %J/s, enthalpy released is dependent on number of moles reacted
    
    %NASA:
    coeff = [5.44417817E+00 -5.20451694E-03 9.62906577E-05 -1.20068814E-07 4.681948255E-11; ...
        5.57278967E+00 3.76541017E-03 6.52226708E-05 -8.30909522E-08 3.20311342E-11; ...
        4.42674073E+00 6.63946249E-03 6.80652815E-05 -9.28753562E-08 3.73473949E-11; ...
        4.19864056E+00 -2.03643410E-03 6.52040211E-06 -5.48797062E-09 1.77197817E-12; ...  
        5.03930607E+00 4.09387100E-04 9.15574112E-05 -1.19411713E-07 4.75043987E-11; ...  
        1.68530424E+00 1.96120012E-02 4.46523571E-05 -8.31523114E-08 3.80651226E-11; ...
        4.45479276E+00 8.26057985E-03 8.29886664E-05 -1.146476425E-07 4.64570101E-11];    
    a1 = coeff(:,1);
    a2 = coeff(:,2);
    a3 = coeff(:,3);
    a4 = coeff(:,4);
    a5 = coeff(:,5);
    Cp = R*(a1+a2*T+a3*T^2+a4*T^3+a5*T^4); %J/mol K 
    
    CpT = Cp'.*x; %J/s K; Cp transpose bc 7x1 * 1x7
    
    % Cooling Duty
    TubeOuterDiameter = TubeInnerDiameter+TubeThickness;
    U = 250; %W/m^2 K, 250-750 for organic solvents cooled by 
    CoolantMassFlow = 9;%kg/s
    CoolantMolarFlow = CoolantMassFlow/0.018; %mol/s
    dHXArea = NumberOfTubes*pi*TubeOuterDiameter*dz;
    dQS = U*dHXArea*(T-CoolantTemp); %J/s 
    %Perry's:
    CpWater = 276.37-2.0901*T+8.125E-03*T^2-0.014116E-03*T^3+9.3701E-09*T^4;
    CoolantTemp = CoolantTemp + dQS/(CpWater*CoolantMolarFlow);
    dT = (-sum(DH)-dQS)/(sum(CpT)*FT);
    %% Pressure drop
    %Ergun
    VolumetricFlow = FT*R*T/(P*1000); %IG Law -- appropriate at high Pressure?
    SuperficialVelocity = VolumetricFlow/FlowArea; 
    Mr = [56.11 56.11 56.11 18.02 74.12 54.09 58.12]./1000; %kg/mol 
    FluidDensity = dot(Mr,x)*P/(R*T);   %kg/m3, assuming IG
    
    %Viscosity coefficients from Riazi pp331-335
    A = [1.0320E-06 1.0320E-06 1.0320E-06 6.1842E-07 1.99E-07 2.6963E-07 6.9154E-07];
    B = [4.8960E-01 4.8960E-01 4.8960E-01 6.7780E-01 7.2330E-01 6.7150E-01 5.2140E-01];
    C = [3.4739E+02 3.4739E+02 3.4739E+02 8.4722E+02 1.7800E+02 1.3472E+02 2.2900E+02];
    D = [0.0000E+00 0.0000E+00 0.0000E+00 -7.4074E+04 0.0000E+00  0.0000E+00 0.0000E+00];
    MuLow = 1000.*A.*T.^B./(1+C./T+D./T^2); %cP; 1cP = 1 mPa s; Pr<0.6 (low pressure)
    %SBA coefficients are actually those for i-propanol
    %SBA data was not available in Riazi, so data from Perry's was used to calculate the vapour viscosity of SBA
    C1 = 1.2114E-07;
    C2 = 0.76972;
    C3 = 92.661;
    MuLow(5) = C1*T^C2/(1+C3/T)*10^3; %cP
    
    Pc = [4043 4043 4043 22055 4302 4277 3648]; %kPa
    Tc = [419.8 419.8 419.8 646.98 547.63 425.02 407.99]; %K
    Pr = P./Pc;
    Tr = T./Tc;
    ReducedDensity = Mr.*Pr./(R.*Tr);
    Zeta = Tc.^(1/6).*Mr.^(-1/2).*(0.987.*Pc).^(-2/3);
   
    Mu = MuLow+((0.1023+0.023364.*ReducedDensity+0.058533.*ReducedDensity.^2-0.040758.*ReducedDensity.^3+0.0093324.*ReducedDensity.^4).^4-10^-4)./Zeta; %cP
    MuMixture = 10^-3*sum(x.*Mr.^0.5.*Mu)/dot(x,Mr.^0.5); %Pa s
    
    dPE = 10^-3*dz*((150*MuMixture*(1-Voidage)^2*SuperficialVelocity)/(Voidage^3*ParticleDiameter^2)+(1.75*(1-Voidage)*FluidDensity*SuperficialVelocity^2)/(Voidage^3*ParticleDiameter)); %kPa 
    dP = -dPE+sum(dF)*R*T/(VolumetricFlow*1000);
    %% Recalculate values
    T = T+dT;
    P = P+dP;
    if P<100
        errordlg('Pressure too low')
        break
    end
    F = F+dF;
    
    FT = sum(F);
    x = F./FT;
    if any(x<0) || any(x>1) 
        errordlg('Unphysical composition')
        break
    end
    X = (sum(F0(1:3))-sum(F(1:3)))/sum(F0(1:3)); %conversion of mixed butenes
    X1 = (F0(1)-F(1))/F0(1);
    X2 = (F0(2)-F(2))/F0(2);
    X3 = (F0(3)-F(3))/F0(3);
    %% Store values in table
    Table(i+1,1) = i*dz;
    Table(i+1,2) = T;
    Table(i+1,3) = P;
    Table(i+1,4) = x(1);
    Table(i+1,5) = x(2);
    Table(i+1,6) = x(3);
    Table(i+1,7) = x(4);
    Table(i+1,8) = x(5);
    Table(i+1,9) = X1;
    Table(i+1,10) = X2;
    Table(i+1,11) = X3;
    Table(i+1,12) = -sum(DH);
    Table(i+1,15) = P;
    Table(i+1,13) = dQS;
    Table(i+1,14) = CoolantTemp;   
end

LogMeanTemp = ((T0-400)-(T-CoolantTemp))/log((T0-400)/(T-CoolantTemp));
HXArea = NumberOfTubes*pi*TubeOuterDiameter*Height;
Q = U*HXArea*LogMeanTemp;
PressureDrop = P0-P;

%% Plot
z = Table(:,1);
x1 = Table(:,4);
x2 = Table(:,5);
x3 = Table(:,6);
xw = Table(:,7);
xa = Table(:,8);
Temp = Table(:,2);
Press = Table(:,3);
X1 = Table(:,9);
X2 = Table(:,10);
X3 = Table(:,11);
Enthalpy = Table(:,12);
Pressure = Table(:,15);
Cooling = Table(:,13);
WaterTemp = Table(:,14);

figure(1)
plot(z,Temp,z,Pressure,z,WaterTemp)
%axis([0,i*dz,300,600])

figure(2)
plot(z,X1,z,X2,z,X3)
legend('cis-butene','trans-butene','1-butene','Location','southeast')
title('Conversion of Butene')
xlabel('Distance along Reactor / m')
ylabel('Conversion')
ax = gca;
ax.Box = 'off';

figure(3)
plot(z,Enthalpy,z,Cooling)

%% Sensitivity Analysis
% plot(z,x1,z,x2,z,x3,z,xa)   %water composition has not been plotted bc is nearly constant at xw = 0.9
% axis([0,i*dz,0,0.1])
% legend('cis-butene','trans-butene','1-butene','SBA','Location','northwest')
% title('Composition')
% xlabel('Distance along Reactor / m')
% ylabel('Composition along reactor')
% ax = gca;
% ax.Box = 'off';
% %saveas(gcf,'Sensitivity7.png')
% 
% x