%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script is used in the data analysis of ASEN 2002 Aero Lab #2.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
clear all; close all; clc;

%% Read/parse data
data = readInput();


%% Calculations
% TODO: Perform calculations with data


%% Error analysis













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

inputFiles = cell(30,1); % Allocate cell array for file names
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

% Replace all unfilled data with NaN - for averaging purposes
avgGroupData(avgGroupData == 0) = NaN;

% Divide data based on section number (e.g. 011, 012, 013)
sec011Data = avgGroupData(:,:,1:10);
sec012Data = avgGroupData(:,:,11:20);
sec013Data = avgGroupData(:,:,21:30);

% Average data across sections
secAllData(:,:,:,1) = sec011Data;
secAllData(:,:,:,2) = sec012Data;
secAllData(:,:,:,3) = sec013Data;
data = mean(secAllData,4,'omitnan');

% Sort data based on incidence angle
data = vertcat(data(:,:,1), data(:,:,2), data(:,:,3));
data = sortrows(data,23);

%% Return extracted data
outputData = data;

end