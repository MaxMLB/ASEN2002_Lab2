%% Script for the second project
tic
clear;
close all;
clc;
data = readtable("ASEN 2002_1230_1_WTData.csv.xlsx");
%waterData = readmatrix("Wind Tunnel Lab Data.csv");
VenturiData = [0.1;0.5;1.4;3;5];

%% Reading in all the data
blfiles = dir('Exp 1 (voltage to speed) data files');
blfiles = blfiles(4:end);
names = cell(20,1); %to hold the names of each file

Data = zeros(100,30,20); %Initialize matrix to store all of the new data

for i= 1:20
    names{i} = strcat('Exp 1 (voltage to speed) data files/',blfiles(i).name);
    Data(:,:,i) = readmatrix(names{i});
end



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



t_avg = averageData(:,1);
p_avg = averageData(:,2);

v_volt = sqrt(2*p_1.*(R_air.*t_avg./p_avg));
v2_volt = sqrt((2.*p_2*R_air.*t_avg)./(p_avg.*(1-(1/9.5)^2)));

%V_volt_Avg = mean([V_volt';V2_volt']);

Voltage1 = averageData(:,13);

figure()
hold on
plot(Voltage1,v_volt);
plot(Voltage1,v2_volt);
hold off


%% Computing the values for each data set and graphing it.

idxVoltage = zeros(20,1);
Voltage = zeros(100,20);
for i = 1:20
    idxVoltage(i) = find(sum(~isnan(Data(:,:,i)),1) > 0, 1 , 'last');
    Voltage(:,i) = Data(:,idxVoltage(i),i);
end
ChangeIndexes = [0;find(diff(Voltage(:,i)) >.1);length(Voltage(:,i))];

AverageData = zeros(length(ChangeIndexes)-1,29,20);

plotVoltage = zeros(5,20);
for i = 1:20
    for j = 1:length(ChangeIndexes)-1
        AverageData(j,:,i) = mean(Data(ChangeIndexes(j)+1:ChangeIndexes(j+1),1:29,i));
        %stdData(i,:) = std(data{ChangeIndexes(i)+1:ChangeIndexes(i+1),1:13});
    end
    plotVoltage(:,i) = AverageData(:,idxVoltage(i),i);
end
P_1 = permute(AverageData(:,5,:),[1,3,2]);
P_2 = permute(AverageData(:,6,:),[1,3,2]);
T_avg = permute(AverageData(:,1,:),[1,3,2]);
P_avg = permute(AverageData(:,2,:),[1,3,2]);



V_volt = sqrt(2*P_1.*(R_air.*T_avg./P_avg));
V2_volt = sqrt((2.*P_2*R_air.*T_avg)./(P_avg.*(1-(1/9.5)^2)));

figure()
plot(plotVoltage,V_volt)
title("Pitostaic Velocity")

figure()
plot(plotVoltage,V2_volt)
title("Venturi Velocity")

toc





