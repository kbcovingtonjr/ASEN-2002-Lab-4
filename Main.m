%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script is used in the data analysis of ASEN 2002 Aero Lab #2.
% -It reads in data taken in a low-speed subsonic wind tunnel including
% pressure distributions across an airfoil.
% -It then calculates axial coefficients and then lift and drag
% coefficients.
% -It then plots this data, similarly to NACA airfoil data.
% 
% Author: Keith Covington
% Created: 11/09/2016
% Modified: 12/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
clear all; close all; clc;


%% Define Airfoil

c = 3.5 * 0.0254; % chord length [m]

% The Pressure Ports coordinates in % Chord
x = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 80, 70, 60, 50, 40, 30, 20,...
    10, 5];
y = [4.19, 9.45, 11.48, 13.6, 14, 13.64, 12.58, 10.95, 8.8, 6.25, 0, 0,...
    0, 0, 0, 0, 0.04, 0.5, 1.11];

skipPositions = [9, 12, 14]; % positions with no pressure data taken

% Remove skipped positions from x and y vectors
x(skipPositions) = [];
y(skipPositions) = [];

% Theoretical port (at trailing edge of airfoil)
theoPortPos = 10;   % position of thoeretical port
theoX = 100;        % x-value at theoretical port
theoY = 0;          % y-value at theoretical port
x = [x(1:theoPortPos-1), theoX, x(theoPortPos:end)]; % add in x-value
y = [y(1:theoPortPos-1), theoY, y(theoPortPos:end)]; % add in y-value

% Scale the Profile for the chord length given
x_scaled = c*x/100;
y_scaled = c*y/100;


%% Read/parse data
inputData = readInput();

% Layer data collected at each speed on top of each other
data(:,:,1) = inputData{1};
data(:,:,2) = inputData{2};
data(:,:,3) = inputData{3};

Patm = data(:,1,:);           % atmospheric pressure [Pa]
Tatm = data(:,2,:);           % atmospheric temperature [K]
density = data(:,3,:);        % atmospheric density [kg/m^3]
speed = data(:,4,:);          % air speed in test section [m/s]
Ppitot = data(:,5,:);         % dynamic pressure [Pa]
Paux = data(:,6,:);           % uhhhhhhhhhh
presPorts = data(:,7:22,:);   % pressures at each port [Pa]
angleDeg = data(:,23,:);      % angle of attack [degrees]
normF = data(:,24,:);         % uhhhhhhhhhh
axialF = data(:,25,:);        % uhhhhhhhhhh
moment = data(:,26,:);        % uhhhhhhhhhh
ELD_X = data(:,27,:);         % uhhhhhhhhhh
ELD_Y= data(:,28,:);          % uhhhhhhhhhh

angleR = angleDeg*pi/180; % changing degrees to radians for calculations

%% Calculations

% Calculate Cp values
Cp1 = findCP(presPorts, Ppitot);

% Extrapolate Cp values to find theoretical Cp value at theoretical port 11
[rows, ~] = size(presPorts);
topLinear = zeros(rows,2);
botLinear = zeros(rows,2);

for speed = 1:length(inputData)
    for angle = 1:rows
        topLinear(angle,:,speed) = polyfit(x_scaled(8:9),...
            Cp1(angle,8:9,speed), 1);
        botLinear(angle,:,speed) = polyfit(x_scaled(10:11),...
            Cp1(angle,10:11,speed), 1);
    end
end

% Evaluate extrapolations at chord length to find theoretical pressures
topEval = topLinear(:,1,:)*c + topLinear(:,2,:);
botEval = botLinear(:,1,:)*c + botLinear(:,2,:);
theoPressures = mean([topEval botEval],2);

% Add theoretical values to Cp matrix
Cp1 = [Cp1(:,1:9,:) theoPressures Cp1(:,10:end,:)];
Cp2 = circshift(Cp1,-1,2);

% Resize some shtuff to work for Cn and Ca calculations
[rows, cols] = size(presPorts);
x1 = repmat(x_scaled,rows,1);
x2 = circshift(x1,-1,2);
y1 = repmat(y_scaled,rows,1);
y2 = circshift(y1,-1,2);

% Cn and Ca calculations
Cn = findCn(Cp1, Cp2, x1, x2, c);
Ca = findCa(Cp1, Cp2, y1, y2, c);

% Cl and Cd calculations
Cl = Cn .* cos(angleR) - Ca .* sin(angleR);
Cd = Cn .* sin(angleR) + Ca .* cos(angleR);

% 0.3205
% 3.002
% 6.9553

% Find maximums on Cl and angle at which it occurs
maxCl = max(Cl)
angleMaxCl = angleDeg(find(Cl == maxCl))

minCd = min(Cd)
angleMinCl = angleDeg(find(Cd == minCd))

%% Plot data

% Plot Cl and Cd
figure

% Plot Cl vs. angle of attack
subplot(1,2,1);
hold on;
plot(angleDeg(:,:,1),Cl(:,:,1),'-*');
plot(angleDeg(:,:,2),Cl(:,:,2),'-*');
plot(angleDeg(:,:,3),Cl(:,:,3),'-*');
title('Lift Coefficient');
xlabel('Angle of Attack (degrees)');
ylabel('Cl');
legend('10 m/s', '20 m/s', '30 m/s', 'location', 'southeast');

% Plot Cd vs. angle of attack
subplot(1,2,2);
hold on
plot(angleDeg(:,:,1),Cd(:,:,1),'-*');
plot(angleDeg(:,:,2),Cd(:,:,2),'-*');
plot(angleDeg(:,:,3),Cd(:,:,3),'-*');
title('Drag Coefficient');
xlabel('Angle of Attack (degrees)');
ylabel('Cd');
legend('10 m/s', '20 m/s', '30 m/s', 'location', 'southeast');


% Plotting Cp graphs
Cp1 = [Cp1 Cp1(:,1,:)];
x_scaled = [x_scaled x_scaled(1)];

angles = [1 9 15 25 30];

for i = 1:length(angles)
    disp(['At ' num2str(angleDeg(angles(i))) ' Degrees:']);
    disp(['Cl is ' num2str(Cl(angles(i),:,:))]);
    disp(['Cd is ' num2str(Cd(angles(i),:,:))]);
    
    figure(i+1)
    hold on;
    plot(x_scaled/c,Cp1(angles(i),:,1),'-o');
    plot(x_scaled/c,Cp1(angles(i),:,2),'-o');
    plot(x_scaled/c,Cp1(angles(i),:,3),'-o');
    set(gca, 'Ydir', 'reverse');
    title(['Cp Over Chord at ' num2str(angleDeg(angles(i))) ' Degrees']);
    legend('10 m/s', '20 m/s', '30 m/s');
    xlabel('Normalized Chord Length');
    ylabel('Cp'); 
end


% Plotting NACA data vs ours
angle_NACA = -8:2:30;
Cl_NACA = [-.13 0.01 0.19 0.33 0.48 0.62 0.78 0.90 1.08 1.19 1.30 1.42...
    1.51 1.60 1.57 1.51 1.45 1.38 1.29 1.14];
Cd_NACA = [0.015 0.014 0.015 0.019 0.022 0.037 0.046 0.059 0.078 0.097...
    0.118 0.14 0.169 0.199 0.242 0.319 0.369 0.416];

figure
subplot(1,2,1);
hold on;
plot(angleDeg(:,:,2),Cl(:,:,2),'-*');
plot(angle_NACA,Cl_NACA,'-*');
title('Lift Coefficient');
xlabel('Angle of Attack (degrees)');
ylabel('Cl');
legend('Our data @ 20 m/s', 'NACA data', 'location', 'southeast');

% Plot Cd vs. angle of attack
subplot(1,2,2);
hold on
plot(angleDeg(:,:,2),Cd(:,:,2),'-*');
plot(angle_NACA(1:end-2),Cd_NACA,'-*');
title('Drag Coefficient');
xlabel('Angle of Attack (degrees)');
ylabel('Cd');
legend('Our data @ 20 m/s', 'NACA data', 'location', 'southeast');



%% Equation for Cn
function Cn = findCn(cp1, cp2, x1, x2, c)
    Cn = 0.5 * (cp1 + cp2) .* (x2 - x1) ./ c;
    Cn = -sum(Cn,2);
end

%% Equation for Ca
function Ca = findCa(cp1, cp2, y1, y2, c)
    Ca = 0.5 * (cp1 + cp2) .* (y2 - y1) ./ c;
    Ca = sum(Ca,2);
end

%% Equation for CP
function Cp = findCP(P, q)
    Cp = P ./ q;
end
