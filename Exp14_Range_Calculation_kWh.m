% ******************************
Velocity=xlsread('FUDS');
V=Velocity(:,2)*1.6;
N=length(V); % No. of Samples
mass = [180 600 1500]; % 2W, 3W & 4W 
area = [0.6 1.6 2.3]; % Frontal area in square metres
Cd = [1 0.45 0.32] ; % Drag coefficient
Gratio = [2/0.28 5/0.2 8/0.3]; % Gearing ratio, = G/r
Pac=[50 200 500]; % Average power of accessories in Watt
Batt_V=[48 48 72];% in Voltage
Batt_Cap=[40 100 200]; %in Ahr
E_Battery=Batt_V.*Batt_Cap/1000;

G_eff = 0.95; % Transmission efficiency
eff_mot=0.80;%Motor Efficiency
Regen_ratio = 1; %This sets the proportion of the

%Initialization of Variables
D_end = zeros(1,100);
Pbat=zeros(1,N);

Pregen=zeros(1,N);
D=zeros(1,N); % Record of distance travelled in km.
Ebat=zeros(1,N);
PTE=zeros(1,N);

i_type=2;% 1 for 2W, 2 for 3W and 3 for 4W
Ebatt_kWh= E_Battery(i_type);%initialize Battery 
SOC_final=100;%initialize (Full Charge Battery)
Distance_total=0;%initial Vehicle stand still

CY=1; %Discharge cycle

while SOC_final>20
for C=2:N
    accel=(V(C) - V(C-1))/3.6;%Acceleration
    D(C) = D(C-1) + (V(C)/3.6/1000);%Distance
    
    Fad = 0.5 * 1.25 * area(i_type) * Cd(i_type) * (V(C)/3.6)^2;% Aero Drag
    Frr=0.015 * mass(i_type) * 9.8; % Rolling Resistance
    Fhc = 0; % Hill Climb Resistance
    Fla = mass(i_type) * accel; % Linear Acceleration
    Pte = (Frr + Fad + Fhc + Fla)*(V(C)/3.6);%Traction Power
    PTE(C)=Pte; % Traction Power
    omega = Gratio(i_type) * (V(C)/3.6);%Motor rotational speed
    
    if omega == 0 % Stationary
        Pte=0;
        Pmot_in=0; % No power into motor
        Torque=0; % No Torque
        
    elseif omega > 0 % Moving        
        if Pte>=0 % Crusing or Accelerating
            Pmot_out=Pte/G_eff; % Motor power > shaft power
            Pmot_in = Pmot_out/eff_mot; % Motoring Mode
        elseif Pte<0 % Braking 
            Pmot_out= Regen_ratio * Pte * G_eff; % Motor power < shaft power
            Pmot_in = Pmot_out * eff_mot; % Generating Mode
        end        
    Torque=Pmot_out/omega; % Torque: Power=Torque*omega
    end
    
    Pbat(C) = Pmot_in + Pac(i_type);   % Power required from Battery
    Ebat(C)=Ebat(C-1)+Pbat(C); % Energy Required from Battery
    
    CY=CY+1;
    ENERGY_B(C)=Ebat(C)/(3600*1000);% Energy consumed so far in this cycle
    ENERGY_C(C)=Ebatt_kWh-Ebat(C)/(3600*1000); %Energy left with battery
    %SOC is Present Charge Available/Total Charge Capacity
    SOC(C)=(Ebatt_kWh-ENERGY_B(C))/Ebatt_kWh*100; % Ignoring Voltage effect

    if SOC_final<=0 % When SOC falls below 0%, then comeout from loop
        SOC_final=0;
        break;
    end
    SOC_final=SOC(C);
%     fprintf('Internal SOC =%f\n',SOC_final);
end
Ebatt_kWh= ENERGY_C(end);

fprintf('SOC at end of each cycle =%f\n',SOC_final);
Distance_total = Distance_total+max(D);
fprintf('Distance travelled at end of each cycle =%f\n',Distance_total);
D=zeros(1,N);
hold on
%Plot for Distance travelled in each cycle vs. SOC
plot(Distance_total,SOC_final,'-*')
xlabel('Distance')
ylabel('SOC')
end
str=['Distance travelled :',num2str(Distance_total),' km'];
title(str);