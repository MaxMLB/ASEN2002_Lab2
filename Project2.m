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

%% Read in the Manometer readings
mano_data = readtable('Water manometer readings.csv');
mano_data_vent = [];
mano_data_pitot = [];

for i = 1:30
   
    if string(mano_data{i,3}) == 'Venturi tube'
        
        mano_data_vent(i,1:10) = [mano_data{i,4:13}];
        mano_data_vent(~any(mano_data_vent,2), :) = [];
        
    elseif string(mano_data{i,3}) == 'Pitot Static Probe'
        
        mano_data_pitot(i,1:10) = [mano_data{i,4:13}];
        mano_data_pitot(~any(mano_data_pitot,2), :) = [];
        
    end
    
end

% get data into voltage in column 1 and height diff in column 2
% more efficient way to do this but im lazy
mano_data_pitot = [mano_data_pitot(:,1:2); mano_data_pitot(:,3:4); mano_data_pitot(:,5:6); mano_data_pitot(:,7:8); mano_data_pitot(:,9:10);];
mano_data_vent = [mano_data_vent(:,1:2); mano_data_vent(:,3:4); mano_data_vent(:,5:6); mano_data_vent(:,7:8); mano_data_vent(:,9:10);];

% calculate differential pressure
% Calculate pressure differential given from manometer
% in. H20 * spec. weight of water [slug/ft^3] * 1/12 [ft] * 32.2 [lb/slug]
% * 47.88 Pa/psf
p_diff_mano_pitot = [mano_data_pitot(:,1), mano_data_pitot(:,2) * 1.940 * (1/12) * 32.2 * 47.88]; % [Pa]
p_diff_mano_vent = [mano_data_vent(:,1), mano_data_vent(:,2) * 1.940 * (1/12) * 32.2 * 47.88]; % [Pa]

% each of the above matrices has the voltage in column 1 and the resultant
% differential pressure in the 2nd column

% Test Plots
%{
subplot(2,1,1);
hold on;
plot(p_diff_mano_pitot(:,1),p_diff_mano_pitot(:,2),'o')
title('Water Manometer - Pitot Tube')
xlabel('Voltage')
ylabel('Pressure Diff [Pa]')
hold off;

subplot(2,1,2);
hold on;
plot(p_diff_mano_vent(:,1),p_diff_mano_vent(:,2),'o')
title('Water Manometer - Venturi Tube')
xlabel('Voltage')
ylabel('Pressure Diff [Pa]')
hold off;
%}

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





