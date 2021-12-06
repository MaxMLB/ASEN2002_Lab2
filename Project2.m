%% Script for the second project
tic
clear;
close all;
clc;
%data = readtable("ASEN 2002_1230_1_WTData.csv.xlsx");
%waterData = readmatrix("Wind Tunnel Lab Data.csv");
R_air = 287; % Gas Constant
%% Reading in all the data
blfiles = dir('Exp 1 (voltage to speed) data files');
blfiles = blfiles(3:end);
names = cell(20,1); %to hold the names of each file

Data = zeros(100,30,20); %Initialize matrix to store all of the new data

for i= 1:20
    names{i} = strcat('Exp 1 (voltage to speed) data files/',blfiles(i).name);
    Data(:,:,i) = readmatrix(names{i});
end

blfiles2 = cell([1,11]);
%load in the data file names
for i = 1:11
    blfiles2{i} = dir(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/BoundaryLayerData/Port ',num2str(i),'/*.csv']);
    %access this by blfiles{1}(1).name for port 1 first file
end

num_elements = zeros(11,1);
for j =1:11
    x=blfiles2{j};
    num_elements(j) = sum(arrayfun(@(x) ~isempty(x.name),x)); %number of datafiles for each port
end
clear i j

blfiles3 = dir('Exp 2 (airfoil) data files');
blfiles3 = blfiles3(2:end);
names3 = cell(25,1); %to hold the names of each file

Data3 = zeros(240,30,25); %Initialize matrix to store all of the new data
blfiles3(1) = [];
for i= 1:25
    names3{i} = strcat('Exp 2 (airfoil) data files/',blfiles3(i).name);
    Data3(:,:,i) = readmatrix(names3{i});
end
%VenturiData = [0.1;0.5;1.4;3;5];
mano_data = readtable('Water manometer readings.csv');

%% Computing the values for each data set

idxVoltage = zeros(20,1);
Voltage = zeros(100,20);
for i = 1:20
    idxVoltage(i) = find(sum(~isnan(Data(:,:,i)),1) > 0, 1 , 'last');
    Voltage(:,i) = Data(:,idxVoltage(i),i);
end
ChangeIndexes = [0;find(diff(Voltage(:,i)) >.1);length(Voltage(:,i))];

AverageData = zeros(length(ChangeIndexes)-1,29,20);
AverageStd = zeros(length(ChangeIndexes)-1,29,20);

plotVoltage = zeros(5,20);
for i = 1:20
    for j = 1:length(ChangeIndexes)-1
        AverageData(j,:,i) = mean(Data(ChangeIndexes(j)+1:ChangeIndexes(j+1),1:29,i));
        AverageStd(j,:,i) = std(Data(ChangeIndexes(j)+1:ChangeIndexes(j+1),1:29,i));
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
mano_data_vent = zeros(50,1);
mano_data_pitot = zeros(50,1);
%formatting and adjusting the csv file matrix
for i = 1:30
    if string(mano_data{i,3}) == "Venturi tube"
        mano_data_vent(i,1:10) = [mano_data{i,4:13}];
        mano_data_vent(~any(mano_data_vent,2), :) = [];
    elseif string(mano_data{i,3}) == "Pitot Static Probe"
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
voltageVals = (0.5:0.5:10);
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
% in. H20 * spec. weight of fluid [slug/ft^3] * 1/12 [ft] * 32.2 [lb/slug]
% * 47.88 Pa/psf
% error in each reading is +/- 0.1 in.
p_diff_mano_pitot = [avg_Mano_Pitot(1,:); avg_Mano_Pitot(2,:) * 1.95 * (1/12) * 32.2 * 47.88]; % [Pa]
p_diff_mano_vent = [avg_Mano_Vent(1,:); avg_Mano_Vent(2,:) * 1.95 * (1/12) * 32.2 * 47.88]; % [Pa]


% each of the above matrices has the voltage in column 1 and the resultant
% differential pressure in the 2nd column

v_mano_pitot = sqrt(2 .* p_diff_mano_pitot(2,:) .* (R_air.*T_avg./P_avg));
v_mano_vent = sqrt((2 .* p_diff_mano_vent(2,:) * R_air .* T_avg)./(P_avg.*(1-(1/9.5)^2)));

volt_diff_mano_pitot = [p_diff_mano_pitot(1,:); p_diff_mano_pitot(1,:); p_diff_mano_pitot(1,:); p_diff_mano_pitot(1,:); p_diff_mano_pitot(1,:)];
volt_diff_mano_vent = [p_diff_mano_vent(1,:);p_diff_mano_vent(1,:);p_diff_mano_vent(1,:);p_diff_mano_vent(1,:);p_diff_mano_vent(1,:)];

[x_hatManoPito] = LSR(v_mano_pitot,volt_diff_mano_pitot);
[x_hatManoVent] = LSR(v_mano_vent,volt_diff_mano_vent);
[sigma_yVentMano,sigma_yPitoMano] = VelocityManometerError(AverageData,plotVoltage,p_diff_mano_pitot',p_diff_mano_vent');

syms X; 
ybestManoPito(X) = x_hatManoPito(1)*X +x_hatManoPito(2);
ybestManoVent(X) = x_hatManoVent(1)*X +x_hatManoVent(2);

figure();
subplot(2,1,1);
hold on;
plot(p_diff_mano_pitot(1,:),v_mano_pitot,'o')
fplot(ybestManoPito,[0 10],"k-","LineWidth",2)
title('Water Manometer - Pitot Tube')
xlabel('Voltage')
ylabel('Velocity [m/s]')
hold off;

subplot(2,1,2);
hold on;
plot(p_diff_mano_vent(1,:),v_mano_vent,'o')
fplot(ybestManoVent,[0 10],"k-","LineWidth",2)
title('Water Manometer - Venturi Tube')
xlabel('Voltage')
ylabel('Velocity [m/s]')
hold off;

%% finding least squares and associated error
Sigma_diffP = 0.01* 6.89476 * 10^3; %Systematic Error in differential pressure transducer
sigmaAbsP = (250-20)*10^3*0.015; %Systematic error in the absoluted pressure values
sigmaT = 0.25; %Systematic error in the temperature data
SigmaPman = 0.05; %Manometer Systematic error

[x_hatPito] = LSR(V_volt,plotVoltage);
[x_hatVent] = LSR(V2_volt,plotVoltage);

[sigmaExtVent,sigmaExtPito] = VelocityError(AverageData,plotVoltage);
yBestPitoVals = (x_hatPito(1).*voltageVals +x_hatPito(2))';
yBestVentVals = (x_hatVent(1).*voltageVals +x_hatVent(2))';
yBestManoVentVals = (x_hatManoVent(1).*voltageVals +x_hatManoVent(2))';

syms X; 
ybestPito(X) = x_hatPito(1)*X +x_hatPito(2);
ybestVent(X) = x_hatVent(1)*X +x_hatVent(2);
figure()
hold on
plot(plotVoltage,V_volt)
fplot(ybestPito,[0 10],"k-","LineWidth",2)
%plot(voltageVals,yBestPitoVals-sigmaExtPito(:,2),"r-","LineWidth",1)
%plot(voltageVals,yBestPitoVals+sigmaExtPito(:,2),"r-","LineWidth",1)
errorbar(voltageVals,yBestPitoVals,sigmaExtPito(:,2));
title("Pitostatic Velocity")
ax = gca;
ax.FontSize = 16;
ylim([0,60])
hold off

figure()
hold on
plot(plotVoltage,V2_volt)
fplot(ybestVent,[0 10],"k-","LineWidth",2)
%plot(voltageVals,yBestVentVals-sigmaExtVent(:,2),"r-","LineWidth",1)
%plot(voltageVals,yBestVentVals+sigmaExtVent(:,2),"r-","LineWidth",1)
errorbar(voltageVals,yBestVentVals,sigmaExtVent(:,2));
title("Venturi Velocity","FontSize",20)
legend("Best Fit Line","Location","northwest","FontSize",20)
xlabel("Voltage [V]","FontSize",18)
ylabel("Velocity [m/s]","FontSize",18)
ax = gca;
ax.FontSize = 16;
ylim([0,60])
hold off


figure()
hold on
fplot(ybestVent,[0 10],"LineWidth",2)
fplot(ybestPito,[0 10],"LineWidth",2)
fplot(ybestManoPito,[0 10],"LineWidth",2)
errorbar(voltageVals,yBestVentVals,sigmaExtVent(:,2));
errorbar(voltageVals,yBestPitoVals,sigmaExtPito(:,2));
errorbar(voltageVals,yBestManoVentVals,sigma_yVentMano(:,2));
title("Different Velocity Measurement Devices","FontSize",20)
legend("Venturi","Pitostatic","Manometer","Venturi Error","PitoStatic Error","Manometer Error","Location","northwest","FontSize",20)
xlabel("Voltage [V]","FontSize",18)
ylabel("Velocity [m/s]","FontSize",18)
ylim([0,60])
ax = gca;
ax.FontSize = 16;

%% Finding the Boundary Layer Thickness
port = struct('meanfreestream',0,'meandensity',0,'blthickness',0,'yloca',zeros(24,1),'speeds',zeros(24,1));
for i = 1:11 %port incrementer
    
    for j = 1:num_elements(i) %number of files in a port incrementer
        
        data = load([blfiles2{i}(j).folder,'/',blfiles2{i}(j).name]);
        %need to group based on data(:,6)
        y = data(:,6);
        y_ind = find(diff(y)>.2); %indicies that correspond to the change of position

        y_ind2 = [0;y_ind;length(y)];
        %and now group the data
        R = 287;
        bobby = NaN(length(y_ind2)-1,8); %preallocate for speeeed!
        for b = 1:length(y_ind2)-1
           temp = mean(data(y_ind2(b)+1:y_ind2(b+1),:)); %the grouped data based on vertical position and then averaged across all samples
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



%% Finding Pressure at trailing edge
angleOfAttack = zeros(4,25);
Pressure = zeros(12,30,25);

for i = 1:25
    for j = 1:12
        Beginidx = (j-1)*20+1;
        EndIdx = j*20;
        Pressure(j,:,i) = mean(Data3(Beginidx:EndIdx,:,i));
        if rem(j,3) == 0
            angleOfAttack(j/3,i) = Pressure(j,7,i);
        end
    end
end
Pressure(:,30,:) = [];

% Finding the average voltage for voltage level 1
VoltVals = zeros(100,3);
x = 1;
for i = 1:25
    for j = 1:4
        VoltVals(x,1) = Pressure((j-1)*3+1,13,i);
        VoltVals(x,2) = Pressure((j-1)*3+2,13,i);
        VoltVals(x,3) = Pressure((j-1)*3+3,13,i);
        x = x+1;
    end
end
clear x
VoltLevels = mean(VoltVals);

Indexes = zeros(25,8);

Indexes(:,1) = angleOfAttack(1,:) == -10;
Indexes(:,2) = angleOfAttack(1,:) == -15;
Indexes(:,3) = angleOfAttack(1,:) == -8;
Indexes(:,4) = angleOfAttack(1,:) == -9;
Indexes(:,5) = angleOfAttack(1,:) == -11;
Indexes(:,6) = angleOfAttack(1,:) == -12;
Indexes(:,7) = angleOfAttack(1,:) == -13;
Indexes(:,8) = angleOfAttack(1,:) == -14;
Indexes = logical(Indexes);
AnglesOfAttack = reshape(angleOfAttack(1:32),[],1); % Grabs all angles of attacks and their occurances
PortPoints = zeros(96,18);
PortPoints(:,1) = [AnglesOfAttack;AnglesOfAttack;AnglesOfAttack];
PortPoints(1:32,2) = VoltLevels(1);
PortPoints(33:64,2) = VoltLevels(2);
PortPoints(65:96,2) = VoltLevels(3);

%Looping through to find a single value for the pressure ports for each
%angle of attack
for i = 1:8
    PortData = Pressure(:,14:29,Indexes(:,i));
    PortDataMean = mean(PortData,3);
    PortPoints((i-1)*4+1:i*4,3:end) = PortDataMean([1,4,7,10],:);
    PortPoints((i-1)*4+33:i*4+32,3:end) = PortDataMean([2,5,8,11],:);
    PortPoints((i-1)*4+65:i*4+64,3:end) = PortDataMean([3,6,9,12],:);
end
clear PortData PortDataMean
PortPoints = sortrows(PortPoints);
% Finding the pressure at the trailing edge from the high and low pressure
% values
P11 = zeros(96,2);
P11(:,1) = PortPoints(:,1);

for j = 1:3
    for i = 1:32
        PHigh = (PortPoints((i-1)*3+j,10)-PortPoints((i-1)*3+j,9)) + PortPoints((i-1)*3+j,10); %Linearly interpolates for the pressure ports on top of the wing
        PLow = (PortPoints((i-1)*3+j,11)-PortPoints((i-1)*3+j,12)) + PortPoints((i-1)*3+j,11); %Linearly interpolates for the pressure ports on bottom of the wing
        P11((i-1)*3+j,2) = mean([PLow,PHigh]); %Finds the average between the 2 linear interpolations and stores it to p11
    end
end
clear PHigh PLow
figure()
P11Graph = P11(1:3:end,:);
plot(P11Graph(:,1),P11Graph(:,2),".-")
title("Pressure at trailing edge for low Velocity")
xlabel("Angle of Attack")
ylabel("Pressure at trailing edge")

figure()
P11Graph = P11(2:3:end,:);
plot(P11Graph(:,1),P11Graph(:,2),".-")
title("Pressure at trailing edge for Medium velocity")
xlabel("Angle of Attack")
ylabel("Pressure at trailing edge")
figure()
P11Graph = P11(3:3:end,:);
plot(P11Graph(:,1),P11Graph(:,2),".-")
title("Pressure at trailing edge for High Velocity")
xlabel("Angle of Attack")
ylabel("Pressure at trailing edge")

%% Calculating Coefficient of Pressure at each point for each angle of attack

P_atm  = mean(Pressure(:,2,:),"all"); %Finds the average freestream pressure
T_atm = mean(Pressure(:,1,:),"all"); %Finds the average freestream temperature
V_free = zeros(3,1);

 %Finds the average freestream velocity from the voltages by average the
 %best fit for the venturi and pito tube
for i = 1:3
    pito = x_hatPito(1)*VoltLevels(i) + x_hatPito(2);
    vent = x_hatVent(1)*VoltLevels(i) + x_hatVent(2);
    V_free(i) =  mean([pito,vent]);
end


rho_free = P_atm ./ (R_air * T_atm); %Calculates freestream density


PortPressure = [PortPoints(:,1:10), P11(:,2), PortPoints(:,11:end)]; %Creates array holding the pressure at each port for each angle of attack for each velocity

C_p = zeros(96,19);
C_p(:,1) = PortPressure(:,1);
C_p(1:3:end,2) = V_free(1);
C_p(2:3:end,2) = V_free(2);
C_p(3:3:end,2) = V_free(3);
for i = 1:96
    if rem(i,3) == 1
        Q_free = 0.5*rho_free*V_free(1)^2; %Calculates freestream dynamic pressure Pa
    elseif rem(i,3) == 2
        Q_free = 0.5*rho_free*V_free(2)^2; %Calculates freestream dynamic pressure Pa
    elseif rem(i,3) == 0
        Q_free = 0.5*rho_free*V_free(3)^2; %Calculates freestream dynamic pressure Pa
    end
    for j = 3:19
        C_p(i,j) = (PortPressure(i,j))/Q_free;
    end
end

%% Plotting the coefficient of pressure vs normalized chord length
xDist = [0;0.175;0.35;0.7;1.05;1.4;1.75;2.1;2.8;3.5;2.8;2.1;1.4;1.05;0.7;0.35;0.175;0]/3.5; %Normalized X chord distanced for each pressure port

%Making a 3d plot of the coeficient of pressure vs angle of attack over the
%entire wing for various velocities
figure()
sgtitle("Coefficient of Pressure over the wing for various angles of attack and velocity")
C_pGraph = [C_p(1:3:end,3:19),C_p(1:3:end,3)];
C_pGraph2 = flip(C_pGraph(:,10:18),2);
C_pGraphPort15 = mean(C_pGraph2(:,6:7),2);
C_pGraph2 = [C_pGraph2(:,1:6),C_pGraphPort15,C_pGraph2(:,7:end)];
C_pGraph = C_pGraph';
C_pGraph2 = C_pGraph2';

subplot(2,2,1)
surf(C_pGraph(1:10,:),'FaceLighting','gouraud','MeshStyle','column','SpecularColorReflectance',0,'SpecularExponent',5,...
    'SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on','LineWidth',1.5,...
    'FaceAlpha',0.5,'FaceColor',[0.07 0.6 1],'EdgeAlpha',1)
hold on
surf(C_pGraph2,'SpecularExponent',1,'MeshStyle','column','SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,...
    'FaceColor',[0.5 0.5 .5],'AlignVertexCenters','on','LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',1)
ylim([1,10])
yticks((1:10))
yticklabels(xDist(1:10))
title("Low Velocity")
ylabel("Normalized Chord Length")
xlim([1,32])
xticks(1:32)
xticklabels(sortrows(AnglesOfAttack))
xlabel("Angle of Attack")
zlabel("Pressure Coefficient")
ax = gca;
ax.ZDir = 'reverse';
ax.Interactions = rotateInteraction;
legend("Top of Wing","Bottom of Wing","Location","northwest")
hold off


C_pGraph3 = [C_p(2:3:end,3:19),C_p(2:3:end,3)];
C_pGraph4 = flip(C_pGraph3(:,10:18),2);
C_pGraphPort15_2 = mean(C_pGraph4(:,6:7),2);
C_pGraph4 = [C_pGraph4(:,1:6),C_pGraphPort15_2,C_pGraph4(:,7:end)];
C_pGraph3 = C_pGraph3';
C_pGraph4 = C_pGraph4';
subplot(2,2,2)
surf(C_pGraph3(1:10,:),'FaceLighting','gouraud','MeshStyle','column','SpecularColorReflectance',0,'SpecularExponent',5,...
    'SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on','LineWidth',1.5,...
    'FaceAlpha',0.5,'FaceColor',[0.07 0.6 1],'EdgeAlpha',1)
hold on
surf(C_pGraph4,'SpecularExponent',1,'MeshStyle','column','SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,...
    'FaceColor',[0.5 0.5 .5],'AlignVertexCenters','on','LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',1)
ylim([1,10])
yticks((1:10))
yticklabels(xDist(1:10))
title("Medium Velocity")
ylabel("Normalized Chord Length")
xlim([1,32])
xticks(1:32)
xticklabels(sortrows(AnglesOfAttack))
xlabel("Angle of Attack")
zlabel("Pressure Coefficient")
ax = gca;
ax.ZDir = 'reverse';
ax.Interactions = rotateInteraction;
legend("Top of Wing","Bottom of Wing","Location","northwest")
hold off

C_pGraph5 = [C_p(3:3:end,3:19),C_p(3:3:end,3)];
C_pGraph6 = flip(C_pGraph5(:,10:18),2);
C_pGraphPort15_3 = mean(C_pGraph6(:,6:7),2);
C_pGraph4 = [C_pGraph6(:,1:6),C_pGraphPort15_3,C_pGraph6(:,7:end)];
C_pGraph5 = C_pGraph5';
C_pGraph6 = C_pGraph6';
subplot(2,2,3)
surf(C_pGraph5(1:10,:),'FaceLighting','gouraud','MeshStyle','column','SpecularColorReflectance',0,'SpecularExponent',5,...
    'SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on','LineWidth',1.5,...
    'FaceAlpha',0.5,'FaceColor',[0.07 0.6 1],'EdgeAlpha',1)
hold on
surf(C_pGraph6,'SpecularExponent',1,'MeshStyle','column','SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,...
    'FaceColor',[0.5 0.5 .5],'AlignVertexCenters','on','LineWidth',1.5,'FaceAlpha',0.5,'EdgeAlpha',1)
ylim([1,10])
yticks((1:10))
yticklabels(xDist(1:10))
title("High Velocity")
ylabel("Normalized Chord Length")
xlim([1,32])
xticks(1:32)
xticklabels(sortrows(AnglesOfAttack))
xlabel("Angle of Attack")
zlabel("Pressure Coefficient")
ax = gca;
ax.ZDir = 'reverse';
ax.Interactions = rotateInteraction;
legend("Top of Wing","Bottom of Wing","Location","northwest")
hold off


subplot(2,2,4)
surf(C_pGraph(1:10,:),'FaceLighting','gouraud','MeshStyle','column','SpecularColorReflectance',0,'SpecularExponent',5,...
    'SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on','LineWidth',1.5,...
    'FaceAlpha',0.8,'FaceColor',[0 0.25 1],'EdgeAlpha',1)
hold on
surf(C_pGraph3(1:10,:),'FaceLighting','gouraud','MeshStyle','column','SpecularColorReflectance',0,'SpecularExponent',5,...
    'SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on','LineWidth',1.5,...
    'FaceAlpha',0.6,'FaceColor',[0 0.863 0.259],'EdgeAlpha',1)
surf(C_pGraph5(1:10,:),'FaceLighting','gouraud','MeshStyle','column','SpecularColorReflectance',0,'SpecularExponent',5,...
    'SpecularStrength',1,'DiffuseStrength',1,'AmbientStrength',0.4,'AlignVertexCenters','on','LineWidth',1.5,...
    'FaceAlpha',0.3,'FaceColor',[1 0.306 0],'EdgeAlpha',1)

ylim([1,10])
yticks((1:10))
yticklabels(xDist(1:10))
title("Pressure over the top for different velocities")
ylabel("Normalized Chord Length")
xlim([1,32])
xticks(1:32)
xticklabels(sortrows(AnglesOfAttack))
xlabel("Angle of Attack")
zlabel("Pressure Coefficient")
ax = gca;
ax.ZDir = 'reverse';
ax.Interactions = rotateInteraction;
legend("Low Velocity","Medium Velocity","High Velocity","Location","northwest")
hold off



% Plotting 2d graph of pressure coefficient over the wing for an angle of
% attack of 3 degrees
C_pGraph2D = [C_p(55:57,3:19),C_p(55:57,3)]';
figure()
plot(xDist,C_pGraph2D,"LineWidth",2)
ax = gca;
ax.YDir = 'reverse';
legend(string(V_free + " m/s"),"Location","northeast","FontSize",20)
title("Coefficient of Pressure over the wing","FontSize",20)
xlabel("Normalized chord length","FontSize",18)
ylabel("Pressure Coefficient","FontSize",18)
ax = gca;
ax.FontSize = 16; 
%% Calculating Normal and axial force coefficient

Cn = zeros(96,1); %Initializing vector
Ca = zeros(96,1);

yDist = [.14665,.33075,.4018,.476,.49,.4774,.4403,.38325,.21875,0,0,0,0,0,.0014,.0175,.03885,0.14665]'/3.5; %Normalized Y location for each pressure port

for i = 1:96
    C_n_tot = zeros(17,1); %Holds the value to be summed for each row
    C_a_tot = zeros(17,1);
    for j = 3:19
        if j ~= 19
            Cpi1 = C_p(i,j+1);
        else
            Cpi1 = C_p(i,3);
        end
        Cpi = C_p(i,j);
        k = j-2;
        DeltaX = xDist(k+1)-xDist(k);
        DeltaY = yDist(k+1)-yDist(k);
        CnElement = 0.5*(Cpi+Cpi1)*DeltaX;
        CaElement = 0.5*(Cpi+Cpi1)*DeltaY;
        C_n_tot(k) = CnElement;
        C_a_tot(k) = CaElement;
    end
    CnLocal = -1* sum(C_n_tot);
    Cn(i) = CnLocal;

    CaLocal = sum(C_a_tot);
    Ca(i) = CaLocal;
end


%% Finding Coefficient of lift and Drag

Cl = Cn.*cosd(C_p(:,1)) - Ca.*sind(C_p(:,1)); %Coefficient of Lift for each angle of attack and velocity
Cd = Cn.*sind(C_p(:,1)) + Ca.*cosd(C_p(:,1)); %Coefficient of drag for each angle of attack and velocity

%Plotting the coefficient of lift for varying velocities and angles of
%attack
figure()
plot(sortrows(AnglesOfAttack),Cl(1:3:end),"LineWidth",1.5)
hold on
plot(sortrows(AnglesOfAttack),Cl(2:3:end),"LineWidth",1.5)
plot(sortrows(AnglesOfAttack),Cl(3:3:end),"LineWidth",1.5)
title("Coefficient of Lift for Different Freestream Velocities","FontSize",20)
xlabel("Angle of Attack [°]","FontSize",18)
ylabel("Coefficient of Lift","FontSize",18)
legend(string(V_free + " m/s"),"FontSize",20,"Location","northwest")
xlim([-15,16])
ax = gca;
ax.FontSize = 16; 
hold off

%Plotting the coefficient of drag for varying velocities and angles of
%attack
figure()
plot(sortrows(AnglesOfAttack),Cd(1:3:end),"LineWidth",1.5)
hold on
plot(sortrows(AnglesOfAttack),Cd(2:3:end),"LineWidth",1.5)
plot(sortrows(AnglesOfAttack),Cd(3:3:end),"LineWidth",1.5)
title("Coefficient of Drag for Different Freestream Velocities","FontSize",20)
xlabel("Angle of Attack [°]","FontSize",18)
ylabel("Coefficient of Drag","FontSize",18)
legend(string(V_free + " m/s"),"FontSize",20,"Location","northwest")
xlim([-15,16])
ax = gca;
ax.FontSize = 16; 
hold off



%% Functions
toc

% Calculating error for the automated sensor data
% Takes in the average data matrix and returns the appropriate least
% squares lines and error
function [sigma_yVent,sigma_yPito] = VelocityError(Data,voltage)
    R = 287; %universal gas constant
    Sigma_diffP = 0.01* 6.89476 * 10^3; %Systematic Error in differential pressure transducer
    sigmaAbsP = (250-20)*10^3*0.015; %Systematic error in the absoluted pressure values
    sigmaT = 0.25; %Systematic error in the temperature data
    areaRatio = (1/9.5);

    PitoError = zeros(5,20);
    VentError = zeros(5,20);
    sigma_yVent = zeros(20,2);
    sigma_yPito = zeros(20,2);
    VoltOrder = [1;3;5;7;9;2;4;6;8;10;1.5;3.5;5.5;7.5;9.5;0.5;2.5;4.5;6.5;8.5];
    sigma_yPito(:,1) = VoltOrder;
    sigma_yVent(:,1) = VoltOrder;
    for i = 1:20
        for j = 1:5
            T_atm = Data(j,1,i);
            P_atm = Data(j,2,i);
            PT1 = Data(j,5,i);
            PT2 = Data(j,6,i);

            %Computing the Error for the pitostatic probe
            d_deltaP = ((R*T_atm)/P_atm)*(((2*PT1*R*T_atm)/(P_atm))^(-1/2));
            d_P= ((-PT1*R*T_atm)/(P_atm^2))*(((2*PT1*R*T_atm)/(P_atm))^(-1/2));
            d_T = ((PT1*R)/(P_atm))*(((2*PT1*R*T_atm)/(P_atm))^(-1/2));
            PitoError(j,i) = sqrt((d_deltaP*Sigma_diffP)^2+(d_P*sigmaAbsP)^2+(d_T*sigmaT)^2);

            %Computing the Error for the Venturi tube
            d_delta_P = ((R*T_atm)/(P_atm*(1-areaRatio^2)))*(((2*PT2*R*T_atm)/(P_atm*(1-areaRatio^2)))^(-1/2));
            dP = ((-PT2*R*T_atm)/((P_atm*(1-areaRatio^2))^2))*(((2*PT2*R*T_atm)/(P_atm*(1-areaRatio^2)))^(-1/2));
            dT = ((PT2*R)/(P_atm*(1-areaRatio^2)))*(((2*PT2*R*T_atm)/(P_atm*(1-areaRatio^2)))^(-1/2));
            VentError(j,i) = sqrt(((d_delta_P * Sigma_diffP)^2) + ((dP * sigmaAbsP)^2) + ((dT * sigmaT)^2));
        end
    end

    Indexes = zeros(20,4);
    Indexes(:,1) = voltage(1,:) == 1;
    Indexes(:,2) = voltage(1,:) == 2;
    Indexes(:,3) = voltage(1,:) == 1.5;
    Indexes(:,4) = voltage(1,:) == 0.5;
    Indexes = logical(Indexes);

    for i = 1:4
        PitoData = PitoError(:,Indexes(:,i));
        VentData = VentError(:,Indexes(:,i));
        sigma_yPito((i-1)*5+1:i*5,2) = mean(PitoData,2);
        sigma_yVent((i-1)*5+1:i*5,2) = mean(VentData,2);
    end
    sigma_yPito = sortrows(sigma_yPito);
    sigma_yVent = sortrows(sigma_yVent);
end


function [sigma_yVent,sigma_yPito] = VelocityManometerError(Data,voltage,pitoDiff,VentDiff)
    R = 287; %universal gas constant
    sigmaAbsP = (250-20)*10^3*0.015; %Systematic error in the absoluted pressure values
    sigmaT = 0.25; %Systematic error in the temperature data
    SigmaPman = 0.05; %Manometer Systematic error
    areaRatio = (1/9.5);

    PitoError = zeros(5,20);
    VentError = zeros(5,20);
    sigma_yVent = zeros(20,2);
    sigma_yPito = zeros(20,2);
    VoltOrder = [1;3;5;7;9;2;4;6;8;10;1.5;3.5;5.5;7.5;9.5;0.5;2.5;4.5;6.5;8.5];
    sigma_yPito(:,1) = VoltOrder;
    sigma_yVent(:,1) = VoltOrder;
    for i = 1:20
        for j = 1:5
            
            T_atm = Data(j,1,i);
            P_atm = Data(j,2,i);
            rho_atm = Data(j,3,i);
            volt = Data(j,13,i);
            diffPPito = pitoDiff(find(pitoDiff(:,1) == volt,1,"first"),2);
            diffPVent = VentDiff(find(VentDiff(:,1) == volt,1,"first"),2);

            %Computing the Error for the pitostatic probe
            DeltaP = 0.5 * ((2/rho_atm)^(0.5)) * (diffPPito)^(-1/2);
            dP = -0.5*((2*diffPPito*R*T_atm)^(0.5))*(P_atm)^(-3/2);
            dT = 0.5* (((2*diffPPito*R)/P_atm)^0.5) * (T_atm^(-0.5));
            PitoError(j,i) = sqrt((DeltaP*SigmaPman)^2+(dP*sigmaAbsP)^2+(dT*sigmaT)^2);

            %Computing the Error for the Venturi tube
            firstPt = (2*R*T_atm) / (P_atm*(1-((areaRatio)^2)));
            partialDeltP = (firstPt^0.5)*0.5*(diffPVent^(-0.5));
            
            firstPt = (2*diffPVent*R*T_atm)/(1-((areaRatio)^2));
            partialP = (firstPt^0.5)*(-0.5)*(P_atm^(-3/2));
            
            firstPt = (2*diffPVent*R)/(P_atm*(1-((areaRatio)^2)));
            partialT = (firstPt^0.5)*(0.5*(T_atm^(-1/2)));
            
            VentError(j,i) = sqrt((partialDeltP*SigmaPman)^2 + (partialT*sigmaT)^2 + (partialP*sigmaAbsP)^2);


%             Delta_P = ((2*R*T_atm) / (P_atm*(1-((areaRatio)^2)))^0.5)*0.5*(diffPVent^(-0.5));
%             d_P = ((2*diffPVent*R*T_atm)/(1-((areaRatio)^2))^0.5)*(-0.5)*(P_atm^(-3/2));
%             d_T = ((2*diffPVent*R)/(P_atm*(1-((areaRatio)^2)))^0.5)*(0.5*(T_atm^(-1/2)));
%             VentError(j,i) = sqrt(((Delta_P * SigmaPman)^2) + ((d_P * sigmaAbsP)^2) + ((d_T * sigmaT)^2));
        end
    end

    Indexes = zeros(20,4);
    Indexes(:,1) = voltage(1,:) == 1;
    Indexes(:,2) = voltage(1,:) == 2;
    Indexes(:,3) = voltage(1,:) == 1.5;
    Indexes(:,4) = voltage(1,:) == 0.5;
    Indexes = logical(Indexes);

    for i = 1:4
        PitoData = PitoError(:,Indexes(:,i));
        VentData = VentError(:,Indexes(:,i));
        sigma_yPito((i-1)*5+1:i*5,2) = mean(PitoData,2);
        sigma_yVent((i-1)*5+1:i*5,2) = mean(VentData,2);
    end
    sigma_yPito = sortrows(sigma_yPito);
    sigma_yVent = sortrows(sigma_yVent);
end



%Least squares estimation of Pitostatic Velocity
function [x_hat] = LSR(y,t)
    d = reshape(y,[],1);
    Acol1 = reshape(t,[],1);
    ALength = length(Acol1);
    Acol2 = ones(ALength,1);
    A = [Acol1,Acol2];

    %computing the x_hat vector
    x_hat = (A'*A)^(-1)*A'*d;
end
