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

%% Read in the Manometer readings
% I moved this down to access the temperature and pressure vals
mano_data = readtable('Water manometer readings.csv');
mano_data_vent = [];
mano_data_pitot = [];

%formatting and adjusting the csv file matrix
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

%{
% Test values
T_avg = 300;
R_air = 287;
P_avg = 101325;
%}

v_mano_pitot = [mano_data_pitot(:,1), sqrt(2 * p_diff_mano_pitot(:,2) .* (R_air.*T_avg./P_avg))];
v_mano_vent = [mano_data_vent(:,1), sqrt((2 .* p_diff_mano_vent(:,2) * R_air .* T_avg)./(P_avg.*(1-(1/9.5)^2)))];

subplot(2,1,1);
hold on;
plot(p_diff_mano_pitot(:,1),v_mano_pitot(:,2),'o')
title('Water Manometer - Pitot Tube')
xlabel('Voltage')
ylabel('Velocity [m/s]')
hold off;

subplot(2,1,2);
hold on;
plot(p_diff_mano_vent(:,1),v_mano_vent(:,2),'o')
title('Water Manometer - Venturi Tube')
xlabel('Voltage')
ylabel('Velocity [m/s]')
hold off;


%% finding least squares and associated error
Sigma_diffP = 0.01* 6.89476 * 10^3;
sigmaAbsP = (250-20)*10^3*0.015;
sigmaT = 0.25;

%finding the uncertainty of each value

WPito  = zeros(100,1);
x = 1;
for i = 1:20
    for j = 1:5
        uncertainty = GUPito(AverageData(j,5,i),AverageData(j,1,i),AverageData(j,2,i));
        WPito(x) = uncertainty;
        x = x+1;
    end
end




[x_hatPito, sigma_yPito] = LSR(V_volt,plotVoltage);
[x_hatVent, sigma_yVent] = LSR(V2_volt,plotVoltage);

syms X; 
ybestPito(X) = x_hatPito(1)*X +x_hatPito(2);
ybestVent(X) = x_hatVent(1)*X +x_hatVent(2);
figure()
hold on
plot(plotVoltage,V_volt)
fplot(ybestPito,[0 10],"k-","LineWidth",2)
fplot(ybestPito-sigma_yPito,[0 10],"r-","LineWidth",1)
fplot(ybestPito+sigma_yPito,[0 10],"r-","LineWidth",1)
title("Pitostaic Velocity")
hold off

figure()
hold on
plot(plotVoltage,V2_volt)
fplot(ybestVent,[0 10],"k-","LineWidth",2)
fplot(ybestVent-sigma_yVent,[0 10],"r-","LineWidth",1)
fplot(ybestVent+sigma_yVent,[0 10],"r-","LineWidth",1)
title("Venturi Velocity")
hold off


toc
%% Least squares estimation of Pitostatic Velocity
function [x_hat,sigma_y_SSE] = LSR(y,t)
    d = reshape(y,[],1);
    Acol1 = reshape(t,[],1);
    ALength = length(Acol1);
    Acol2 = ones(ALength,1);
    A = [Acol1,Acol2];
    %computing the x_hat vector
    x_hat = (A'*A)^(-1)*A'*d;
    diff_y = d - (x_hat(1)*Acol1 + x_hat(2)); %difference between the actual data and best fit line
    sigma_y_SSE = sqrt((1/(length(diff_y)-2))*sum(diff_y.^2)); %SSE error from the best fit line
end

function result = GUPito(Delta_P,T,P)
    Sigma_diffP = 0.01* 6.89476 * 10^3;
    sigmaAbsP = (250-20)*10^3*0.015;
    sigmaT = 0.25;
    R = 0.287;
    d_deltaP = (2 .* R .* T./P) .* (Delta_P.*R.*T./P).^(-1/2);
    d_T = (2.*Delta_P.*R./P) .* (Delta_P.*R.*T./P).^(-1/2);
    d_P = (-2*Delta_P.^2.*R^2.*T.^2)/(P.^3);

    result = sqrt((d_deltaP*Sigma_diffP)^2+(d_P*sigmaAbsP)^2+(d_T*sigmaT)^2);
end

