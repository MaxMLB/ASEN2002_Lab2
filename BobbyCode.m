clc; clear all; close all;

%% Wind tunnel calibration

%% Wind tunnel boundary layer


%load in the data file names
for i = 1:11
    blfiles{i} = dir(['Aero Lab Windtunnel Calibration/Aero Lab 1 - 2019 Group Data/BoundaryLayerData/Port ',num2str(i),'/*.csv']);
    %access this by blfiles{1}(1).name for port 1 first file
end

for j =1:11
    x=blfiles{j};
    num_elements(j) = sum(arrayfun(@(x) ~isempty(x.name),x)); %number of datafiles for each port
end

clear i j

for i = 1:11 %port incrementer
    for j = 1:num_elements(i) %number of files in a port incrementer

        data = load([blfiles{i}(j).folder,'/',blfiles{i}(j).name]);
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
    
    %% Now to find the thickness at the port
    
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

%% And now compare to theory
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


% P = polyfit(probespeeds,ylocation,2);
% x = linspace(min(probespeeds),max(probespeeds),100);
% fitline = P(1)*x.^2 + P(2)*x + P(3);
%lets try fit instead
% f = fit(probespeeds,ylocation,'exp1');

%figure;plot(f,probespeeds,ylocation)

%% do the uncertainty analysis

%% look at the x_cr distance of a flat plate
