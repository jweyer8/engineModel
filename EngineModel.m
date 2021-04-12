% Engine Project
% Jared, W
% November 11, 2020

%cleaning
clear all
close all
clc

%range of crack angles
step = .5;               %step size of theta n to n+1
theta = (-180:step:180); %[deg] crank angle 

%engine constants 

%define the rpms
rpms = [1000:100:7000];

%Prepopulate arrays
powerAI = zeros(length(rpms), 1);
powerfAI = zeros(length(rpms), 1);
powerAl = zeros(length(rpms), 1);
powerfAl = zeros(length(rpms), 1);
torqueAI = zeros(length(rpms), 1);
torquefAI = zeros(length(rpms), 1);
torqueAl = zeros(length(rpms), 1);
torquefAl = zeros(length(rpms), 1);

%Loop through rpm values
for i = 1:length(rpms)
    %call engine fucntion with one with  
    [QdotAry, PAry, VAry, TgAry, dVAry]=engine(theta,step,0,rpms(i));
    [QdotAryI, PAryI, VAryI, TgAryI, dVAryI]=engine(theta,step,1,rpms(i));
   
    %get the work done by the piston using ideal (no heat loss)
    [power, powerf, torque, torquef] = work(PAryI, VAryI, rpms(i));
    powerAI(i) = power;
    powerfAI(i) = powerf;
    torqueAI(i) = torque;
    torquefAI(i) = torquef;
    
    %get the work done by the piston using ideal (no heat loss)
    [power, powerf, torque, torquef] = work(PAry, VAry, rpms(i));
    powerAl(i) = power;
    powerfAl(i) = powerf;
    torqueAl(i) = torque;
    torquefAl(i) = torquef;
    
    %plotting durring last iteration 7000 rpm
    if(i == length(rpms))
        figure
        tiledlayout(2,2)
        nexttile
        plot(theta, VAry, 'lineWidth', 2);
        xlim([-180 180])
        title('Volume vs Crank Angle');
        xlabel('Crank Angle [deg]');
        ylabel('Piston Volume [m^3]');
        grid on
        nexttile
        plot(theta, QdotAry, 'lineWidth', 2);
        xlim([-180 180])
        title('Qdot vs Crank Angle');
        xlabel('Crank Angle [deg]');
        ylabel('Heat Release [KJ]');
        grid on
        nexttile
        plot(theta, TgAry, 'lineWidth', 2);
        xlim([-180 180])
        title('Gas Temperature vs Crank Angle');
        xlabel('Crank Angle [deg]');
        ylabel('Gas Temperatute [K]');
        grid on
        nexttile
        plot(theta, PAry,'lineWidth', 2);
        xlim([-180 180]);
        title('Pressure vs Crank Angle');
        xlabel('Crank Angle [deg]');
        ylabel('Piston Pressure [KPa]');
        grid on
        
        %Close PV diagram
        %Nececary becuause not all strokes are modeled 
        VAry(end + 1) = VAry(end);
        VAryI(end + 1) = VAryI(end);
        PAry(end + 1) = 101.325;
        PAryI(end + 1) = 101.325;
        
        %plotting PV diagram
        figure
        hold on
        plot(VAry, PAry, 'lineWidth', 2);
        plot(VAryI, PAryI, 'lineWidth', 2);
        title('PV Diagram');
        ylabel('Pressure [KPa]');
        xlabel('Volume [m^3]');
        legend('Actual', 'Ideal','location',...
               'southoutside','Orientation','horizontal');
        grid on
    end
end

%Fit test data from engine
rpmTests = [2500,3000,3500,4000,5000,5500,5800]; %[rpm] rpm test points
prelimT = [3.25,4.55,5.2,5.59,5.2,4.94,4.55];    %[Nm] Prelim data torques
prelimP = prelimT.*2.*pi.*rpmTests/60;           %[KW] Prelim data power
postT = [1.36,2.71,4.07,4.52,4.52,4.52,3.62];    %[Nm] post data torques
postP = postT.*2.*pi.*rpmTests/60;               %[KW] post data power

%plotting power data
figure
hold on
grid on
plot(rpms, powerAI, '-.', 'lineWidth', 3);
plot(rpms, powerfAI,'--', 'lineWidth', 3);
plot(rpms, powerAl,':', 'lineWidth', 3);
plot(rpms, powerfAl,'-', 'lineWidth', 3);
plot(rpmTests, prelimP/1000,'o', 'lineWidth', 3); 
plot(rpmTests, postP/1000,'o', 'lineWidth', 3);
title('Power vs RPM');
ylabel('Power [KW]');
xlabel('N [rpm]');
legend('Ideal','Friction Loss','Heat Loss', 'Heat + Friction',...
       'Preliminary Data', 'Post Process Data',...
       'location', 'southoutside','Orientation','horizontal', 'NumColumns', 3);
hold off

%plotting torque data
figure
hold on
grid on
plot(rpms, torqueAI*1000,'-.', 'lineWidth', 3);
plot(rpms, torquefAI*1000,'--', 'lineWidth', 3);
plot(rpms, torqueAl*1000,':', 'lineWidth', 3);
plot(rpms, torquefAl*1000,'-', 'lineWidth', 3);
plot(rpmTests, prelimT,'o', 'lineWidth', 3); 
plot(rpmTests, postT,'o', 'lineWidth', 3);
title('Torque vs RPM');
ylabel('Torque [Nm]');
xlabel('N [rpm]');
legend('Ideal', 'Friction Loss', 'Heat Loss', 'Heat + Friction',...
       'Preliminary Data', 'Post Process Data',...
       'location', 'southoutside','Orientation','horizontal','NumColumns', 3);
hold off
    

%function for plotting cylinder pressure as a function of crank angle
%Parameters:
%theta = range of crank angles 
%step = length between theta_n and theta_n+1
%ideal = 1 if no heat loss
%N = rpm of engine
function[QdotAry, P, VAry, TgAry, dVAry] = engine(theta, step, ideal, N)
    %constants
    stdP = 101.325;  %[KPa] Pressure at standard conditions
    stdT = 300;      %[K] temperture when mixture enters piston
    a = 5;           %wiebe function parameter
    n = 3;           %wiebe function parameter
    thetad = 40;     %[deg] combustion duration
    gamma = 1.4;     %specific heat ratio assumed constant
    Mmix = 30.24;    %[kg/kmol] mass of fuel mixture
    AFs = 15.05;     % stoich air fuel ratio
    Ru = 8.314;      %[J/mol*K] Universal gas constant
   
    %engine spec constants measured at lab
    r = 8.5;       %compression ratio, found in the manual
    b = .0397;     %[m] bore of pistion
    s = .0383;     %[m] stroke of the pistion
    L = .06985;    %[m] length of connecting rod
    spark = -18;   %[deg] ignition crank angle

    %solved constants
    Vd = pi*b^2*s/4  %[m^3] displacement volume
    Vc = 8.199e-6;   %[m^3] clearence volume
    c = 2*L/s;  
    
    %calculate qin this is depedent of the iherient heating value of the
    %fuel and the mass of the fuel coming into the cylinder per cycle
    rho = stdP/((Ru/Mmix)*stdT)  %[mols]# mol per cycle
    m = rho*(Vd + Vc)            %[kg] mass of mix per cycle
    mfuel = m/AFs                %[kg] mass of fuel per cycle
    HV = 43000;                  %[kJ] iso-octane fuel
    qin = mfuel*HV;              %[kJ/kg] energy of fuel 
    
    %Prealocate pressure, energy release, volume, volume change rate, and
    %gas temperature 
    P = zeros(1, length(theta));
    QdotAry = zeros(1, length(theta));
    VAry = zeros(1, length(theta));
    dVAry = zeros(1, length(theta));
    TgAry = zeros(1, length(theta));

    %intitial condition at BDC
    P(1) = stdP;

    %solve for values depedent on crank angle
    for i = 1:(length(P)) 
         %wiebe function
         xb = 1-exp(-a*((theta(i)-spark)/thetad)^n);
         
         %Volume of cylinder as a function of crank angle
         V=Vd/(r-1)+Vd/2*(c+1-cosd(theta(i))-sqrt(c^2-(sind(theta(i)))^2)); 
         
         %energy release from flame
         %ideal minus loss to wall due to convection
         %if ideal -> no heat loss
         [Qdotloss, Tg] = heatLoss(N, b, V, Vd, Vc, P(i), gamma, s, m, Mmix);
         if(ideal)
             Qdotloss = 0;
         end
         Qdot=((n*a*qin)/thetad)*(1-xb)*((theta(i)-spark)/thetad)^(n-1);
         Qdottot = Qdot - Qdotloss;
         
         %volume change rate as function of crank angle
         Vdot = Vd*sind(theta(i))*(1 + cosd(theta(i))...
                *(c^2 - sind(theta(i))^2)^(-.5))*(pi/360); 

         %heating term is 0 before spark/combustion
         if ~(theta(i) > spark && theta(i) < spark + thetad)
             Qdottot = 0;
         end

         %eulers method
         %Pressure change rate as a function of theta
         Pdot = ((gamma-1)/V)*(Qdottot) - (gamma/V)*P(i)*Vdot;
         %Pressure in cylinder at theta
         P(i+1) = P(i) + step*Pdot;

         %Fill prealocated array with data for specific theta
         QdotAry(i) = Qdottot;
         VAry(i) = V;
         TgAry(i) = Tg;
         dVAry(i) = Vdot;
    end
    P(end) = []; %Need Pressure array to be same size as others
end


%function for finding the heat transfer out of the cylinder walls
%arguments are: rpm, bore, volume, displacment volume, pressure
%ratio of Cp to Cv, stroke, number of moles entering piston per cycle
function [QdotWall, Tg] = heatLoss(N, b, V, Vd, Vc, P, gamma, s, m, Mmix)
    %reference constants
    Ru = 8.314;   %[J/mol*K] Universal gas constant
    Pr = 101.325; %[KPa] Pressure @ IVC 
    Tr = 300;     %[K] Temperature @ IVC 
    Tw = 500;     %[K] Temperature of wall assumed constant
    Vr = Vc + Vd; %[m^3] Volume of cylinder @IVC
    
    %temperatures of the gas and the cylinder wall
    Tg = (P*V)/((Ru/Mmix)*m); %[K] Temperature of gas assume mixture is ideal
 
    %Calculate the heat loss due to covection
    Up = 2*N*s/60;                               %[m/s] mean piston speed
    Pm = Pr*Vd^gamma/V^gamma;                    %[KPa] Motorized pressure 
    U = (2.28*Up+.00324*Tr*(Vd/Vr)*((P-Pm)/Pr)); %[m/s] Gas velocity
    hg = (3.26*P^.8*U^.8*b^-.2*Tg^-.55)/1000;    %[KW/m^2K] HT coeffient
    A = (2*pi*b^2/4)-(4*Vc/b)+(4*V/b);           %[m^2] Piston chamber area
    QdotWall = hg*A*(Tg - Tw)*(1/(N/60))/360;    %[KJ] Rate of heat loss 
end


%funtion for finding the losses do to friction 
%takes rmp as argument 
function fmep = frictionLoss(N)
    fmep = 0.97 + 0.15*(N/1000) + 0.05*(N/1000)^2; %[bar] friction loss
end


%function for finding the work, power, and torque done by the piston
%arguments are pressure & volume arrays & the rpm
function [power, powerf, torque, torquef] =  work(P, V, N)
    Vd = 4.7410e-05;                  %[m^3] displacment volume 
    work = trapz(V, P);               %[KJ] work done by piston
    power = work*(N/60);              %[KW] power without friction losses
    torque = power*60/(2*pi*N);       %[KNm] torque without friction losses
    bmep = power/(Vd*(N/60));         %[KPa] power in pressue
    fmep = frictionLoss(N)*100;       %[KPa] friction losses in pressure
    powerf = (Vd*(N/60))*(bmep-fmep); %[KW] power with friction loss
    torquef = powerf*60/(2*pi*N);     %[KNm] torque with friction losses
end