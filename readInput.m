function outputData = readInput()
% function data = readInput()
% 
% Reads in the data from the wind tunnel experiements
% performed in ASEN 2002 Aero Lab #2.
%
% input:  --------------
%
% output: outputData  -  matrix of data averaged across sections organized
%                        by ascending value of incidence angle
%
% Author: Keith Covington, Nov 9 2016
% Modified: 11/09/2016

addpath('WindTunnelData');

%% Define input files
inputFiles = cell(10,1); % Allocate cell array for file names
count = 1;

% Generate file names of group test data
for group = 1:9
    inputFiles{count} = ['AirfoilPressure_S013_G0' num2str(group)...
        '_LA.csv'];
    count = count + 1;
end
inputFiles{10} = 'AirfoilPressure_S013_G10_LA.csv'; % they suck at consistent naming
% disp(inputFiles); % Display generated file names


%% Extract data from input files

w = waitbar(0,'Reading wind tunnel data files...'); % progress bar
numFiles = numel(inputFiles); % define number of input files for reference
avgGroupData = zeros(9,28,numFiles); % matrix of each group's average data

% Loop through all data files
for file = 1:numFiles
    
    % Extract data from file if it exists
    try
        % Read .csv file and organize data into cell arrays
        inputData = csvread(inputFiles{file},1); % read .csv file
        [~, cols] = size(inputData);
        cellRowDiv = [500 500 500 500 500 500 500 500 500]; % diminsions for dividing up into cell arrays
        outputData = mat2cell(inputData,cellRowDiv,cols); % organized data
        
        % Average data from each DAQ incidence sent
        for i = 1:9
            avgGroupData(i,:,file) = mean(outputData{i});
        end
                
    % If file doesn't exist, output message and continue
    catch
        disp(['File:  ' inputFiles{file}...
            '  does not exist yet. Skipping file...']);
    end
    
    waitbar(file/numFiles); % update progress bar
end
close(w); % close progress bar


data = avgGroupData(:,:,1:10);

% Sort data based on incidence angle
data = vertcat(data(:,:,1), data(:,:,2), data(:,:,3), data(:,:,4),...
    data(:,:,5), data(:,:,6), data(:,:,7), data(:,:,8), data(:,:,9),...
    data(:,:,10));
data = sortrows(data,23);

% Sort data into 3 arrays based on airspeed
data10 = data(1:3:end,:);   % data collected at 10 m/s
data20 = data(2:3:end,:);   % data collected at 20 m/s
data30 = data(3:3:end,:);   % data collected at 30 m/s

%% Return extracted data
outputData = {data10, data20, data30};



end
