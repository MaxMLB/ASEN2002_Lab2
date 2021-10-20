%Script for the second project
clear;
close all;
clc;
data = readtable("ASEN 2002_1230_1_WTData.csv.xlsx");

R_air = 287;

voltage = data.AppliedVoltage_V_;

ChangeIndexes = [0;find(diff(voltage) >.1);100];


averageData = zeros(5,13);
%stdData = zeros(5,13);


for i = 1:length(ChangeIndexes)-1
    averageData(i,:) = mean(data{ChangeIndexes(i)+1:ChangeIndexes(i+1),1:13});
    %stdData(i,:) = std(data{ChangeIndexes(i)+1:ChangeIndexes(i+1),1:13});
end
delta_p = mean([averageData(:,5)';averageData(:,6)'])';
%delta_pStd = std([averageData(:,5)';averageData(:,6)'])';


T_avg = averageData(:,1);
P_avg = averageData(:,2);

V_volt = sqrt(2*delta_p.*(R_air.*T_avg./P_avg));



