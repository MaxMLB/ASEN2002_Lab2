%Script for the second project
clear;
close all;
clc;
data = readtable("ASEN 2002_1230_1_WTData.csv.xlsx");
%waterData = readmatrix("Wind Tunnel Lab Data.csv");
VenturiData = [0.1;0.5;1.4;3;5];

R_air = 287;

voltage = data.AppliedVoltage_V_;

ChangeIndexes = [0;find(diff(voltage) >.1);length(voltage)];


averageData = zeros(5,13);
%stdData = zeros(5,13);


for i = 1:length(ChangeIndexes)-1
    averageData(i,:) = mean(data{ChangeIndexes(i)+1:ChangeIndexes(i+1),1:13});
    %stdData(i,:) = std(data{ChangeIndexes(i)+1:ChangeIndexes(i+1),1:13});
end
p_1 = averageData(:,5);
p_2 = averageData(:,6);



T_avg = averageData(:,1);
P_avg = averageData(:,2);

V_volt = sqrt(2*p_1.*(R_air.*T_avg./P_avg));
V2_volt = sqrt((2.*p_2*R_air.*T_avg)./(P_avg.*(1-(1/9.5)^2)));

%V_volt_Avg = mean([V_volt';V2_volt']);

Voltage = averageData(:,13);

figure()
hold on
plot(Voltage,V_volt);
plot(Voltage,V2_volt);
%plot(Voltage,V_volt_Avg)






