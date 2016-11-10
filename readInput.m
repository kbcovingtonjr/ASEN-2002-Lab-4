function outputData = readInput()
% function data = readInput()
% 
% Reads in the data from the wind tunnel experiements
% performed in ASEN 2002 Aero Lab #2.
%
% input:  --------------
%
% output: data  - 
%
% Author: Keith Covington, Nov 9 2016
% Modified: 11/09/2016

addpath('WindTunnelData');

%% Define input files

% inputFiles = 'SampleData.csv';

inputFiles = cell(30,1);
count = 1;

% Generate file names of group test data
for sec = 1:3
    for group = 1:10
        inputFiles{count} = ['AirfoilPressure_S01' num2str(sec) '_G0'...
            num2str(group) '.csv'];
        count = count + 1;
    end
end
% disp(inputFiles); % Display generated file names


%% Extract data from input files

w = waitbar(0,'Reading wind tunnel data files...');

% Loop through all data files
for file = 1:numel(inputFiles)
    
    % Extract data from file if it exists
    try
        inputData = csvread(inputFiles{file},1);
        [~, cols] = size(inputData);
        cellRowDiv = [500 500 500 500 500 500 500 500 500];
        outputData = mat2cell(inputData,cellRowDiv,cols);
        
        % Average data from each DAQ incidence sent
        avgData = zeros(9,cols);
        for i = 1:9
            avgData(i,:) = mean(outputData{i});
        end
%         disp(avgData(:,23)); % Display angles of attack from each group
        
        % Sort data
        
        
    % If file doesn't exist, output message and continue
    catch
        disp(['File:  ' inputFiles{file}...
            '  does not exist yet. Skipping file...']);
    end
    
    waitbar(file/numel(inputFiles));
end
close(w);



%% Return extracted data
outputData = 0;
end