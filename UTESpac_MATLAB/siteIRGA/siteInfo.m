% site specific script to load site information

% enter the sonic height in input data table is ascending (true) or descending (false) when multiple sonic heights exist
info.ascending = false;

% enter orientation of sonics.  Sonic order: tables sorted alphabetically followed by columns sorted in descending order
info.sonicOrientation = [243, 128, 134, 139]; %UU1

% enter manufacturer of SAT.  1 for Campbell, 0 for RMYoung.  RMYoung v = Campbell u!
info.sonicManufact = [1, 1, 1, 1];  % 3 IRGASONs

% enter orientation of tower relative to sonic head
info.tower = 210; % FM 

% tower elevation
info.siteElevation = 1980; % (m)

% site slope angle
info.angle = 8.2; % mean of 20 m radius buffer

% enter expected table names.  Missing tables will be filled with NaNs to create consistency 
% when multiple output files are concatenated with getData.m
info.tableNames = {'FMDOL_20Hz', 'FMDOL_1min'}; % modified by Diane

% enter table scan frequencies corresponding to tableNames
info.tableScanFrequency = [20, 1/60];  %[Hz]

% enter number of columns in each .csv table.  Note that the number of columns in the output structure will 
% be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a Matlab's 
% single-column serial time.  Also, note that View Pro frequently cuts of column 1 (the year!) of the .csv file. 
info.tableNumberOfColumns = [48, 10]; % modified by Diane

% use local mean reference temperature from slow sensors if exists
info.useTrefHMP = true;
info.avgSlowFreq = 1; % averaging period to match fast to slow
info.shiftsSonHeight = [13.94, 32.18]; % needs to be exist in pairs with below
info.shiftsHMPHeight = [15, 30];

% seperately run sonic at 51.5 m but need to zRef at 4.42 m
info.shiftzRef = false; % only if run a single sonic at 51.5 m
info.zRefLowestSon = 4.42;

% SSITC quality-flagging settings
info.SSITC_subAvgMin = 5;

% Set this to your best estimate of canopy height.
info.canopyHeight = 19.3;

% Start with 0 for simplicity.
% You can later test 0.67*canopyHeight for above-canopy z-d correction.
info.displacementHeight = 0;

% Use canopy-aware sigma_w/u* ITC model for z <= canopyHeight.
info.useCanopyITC = true;