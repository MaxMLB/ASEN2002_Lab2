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
blfiles = blfiles(3:end);
names = cell(20,1); %to hold the names of each file

Data = zeros(100,30,20); %Initialize matrix to store all of the new data

for i= 1:20
    names{i} = strcat('Exp 1 (voltage to speed) data files/',blfiles(i).name);
    Data(:,:,i) = readmatrix(names{i});
end

%load in the data file names
for i = 1:11
    blfiles2{i} = dir(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/BoundaryLayerData/Port ',num2str(i),'/*.csv']);
    %access this by blfiles{1}(1).name for port 1 first file
end

for j =1:11
    x=blfiles2{j};
    num_elements(j) = sum(arrayfun(@(x) ~isempty(x.name),x)); %number of datafiles for each port
end


clear i j




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
mano_data_vent = zeros(50,1);
mano_data_pitot = zeros(50,1);
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

mano_data_pitot = sort(mano_data_pitot);
mano_data_vent = sort(mano_data_vent);

% calculate average values at each voltage point
voltageVals = 0.5:0.5:10;
avg_Mano_Pitot = [voltageVals; zeros(1,20)];
avg_Mano_Vent = [voltageVals; zeros(1,20)];

for i = 1:length(voltageVals)
    
    idx1 = mano_data_pitot == voltageVals(i);
    idx2 = mano_data_vent == voltageVals(i);
    avg_Mano_Pitot(2,i) = mean(mano_data_pitot(idx1(:,1),2));
    avg_Mano_Vent(2,i) = mean(mano_data_vent(idx2(:,1),2));
    
    
end

% calculate differential pressure
% Calculate pressure differential given from manometer
% in. H20 * spec. weight of water [slug/ft^3] * 1/12 [ft] * 32.2 [lb/slug]
% * 47.88 Pa/psf
p_diff_mano_pitot = [avg_Mano_Pitot(1,:); avg_Mano_Pitot(2,:) * 1.940 * (1/12) * 32.2 * 47.88]; % [Pa]
p_diff_mano_vent = [avg_Mano_Vent(1,:); avg_Mano_Vent(2,:) * 1.940 * (1/12) * 32.2 * 47.88]; % [Pa]

% each of the above matrices has the voltage in column 1 and the resultant
% differential pressure in the 2nd column

v_mano_pitot = sqrt(2 .* p_diff_mano_pitot(2,:) .* (R_air.*T_avg./P_avg));
v_mano_vent = sqrt((2 .* p_diff_mano_vent(2,:) * R_air .* T_avg)./(P_avg.*(1-(1/9.5)^2)));

figure();
subplot(2,1,1);
hold on;
plot(p_diff_mano_pitot(1,:),v_mano_pitot,'o')
title('Water Manometer - Pitot Tube')
xlabel('Voltage')
ylabel('Velocity [m/s]')
hold off;

subplot(2,1,2);
hold on;
plot(p_diff_mano_vent(1,:),v_mano_vent,'o')
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


%% Finding the Boundary Layer Thickness

for i = 1:11 %port incrementer
    for j = 1:num_elements(i) %number of files in a port incrementer

        data = load([blfiles2{i}(j).folder,'/',blfiles2{i}(j).name]);
        %need to group based on data(:,6)
        y = data(:,6);
        y_ind = find(diff(y)>.2); %indicies that correspond to the change of position

        y_ind = [0;y_ind;length(y)];
        %and now group the data
        R = 287;
        bobby = NaN(length(y_ind)-1,8); %preallocate for speeeed!
        for b = 1:length(y_ind)-1
           temp = mean(data(y_ind(b)+1:y_ind(b+1),:)); %the grouped data based on vertical position and then averaged across all samples
           bobby(b,:) = [temp,sqrt(2*temp(4)*(R*temp(2)/temp(1)))]; %put in the probe airspeed as part of the array
           %bobby is an array of the collapsed data (one row per vertical location
           %and the probe airspeed in the 8th column. Note: for one data
           %file only
        end

        offset = .5; %wall thickness of the pitot probe
        freestream = bobby(end,8); % freestream velocity
        ylocation = [0;bobby(1:end-1,6) + offset]; %vertical location
        probespeeds = [0;bobby(1:end-1,8)]; %airspeed sensed by probe

        %group things based on the port number
        
        port(i).test(j).freestream = freestream;
        port(i).test(j).ylocation = ylocation;
        port(i).test(j).probespeeds = probespeeds;
        port(i).test(j).density = temp(1)/(R*temp(2));
    end
    
    % Now to find the thickness at the port
    
    temptable = struct2table(port(i).test);
    port(i).meanfreestream = mean(temptable.freestream);
    port(i).meandensity = mean(temptable.density);
    port(i).yloca = [port(i).test(1).ylocation;port(i).test(2).ylocation];
    port(i).speeds = [port(i).test(1).probespeeds;port(i).test(2).probespeeds];
    
    %figure(i);plot(port(i).speeds,port(i).yloca,'.')
    f = fit(probespeeds,ylocation,'exp1'); %fit some exponential to the data
    port(i).blthickness = f(.99*port(i).meanfreestream);
    %  figure;plot(f,probespeeds,ylocation)
    %  hold on; plot(.95*port(i).meanfreestream,port(i).blthickness,'*')
end

allcombine = struct2table(port);

port_real_locations = 0.0254*[9.05;10.03;11.01;11.99;12.97;13.95;14.93;15.91;16.89;17.87;18.85]; %mm

% And now compare to theory
mu = 1.7894e-5; %kg/(m)(s)
densityinf = mean(allcombine.meandensity);
Vinf = mean(allcombine.meanfreestream);

x_for_bl = linspace(min(port_real_locations),max(port_real_locations),100);
Rex = densityinf*Vinf*x_for_bl/mu;

thick_laminar = 5.2*x_for_bl./sqrt(Rex);
thick_turb = (0.37*x_for_bl)./(Rex.^(0.2));
figure;plot(port_real_locations,allcombine.blthickness,'.')
hold on;
plot(x_for_bl,1000*thick_laminar,'-')
plot(x_for_bl,1000*thick_turb,'--')
xlabel('Port location [m]')
ylabel('Boundary Layer thickness [mm]')
legend('Data','Laminar Theory','Turbulent Theory')

toc

%% Functions

%Least squares estimation of Pitostatic Velocity
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