% site specific script to load site information
% ES5 Spring

% enter orientation of sonics.  Sonic order: tables sorted alphabetically followed by columns sorted in ascending order
info.sonicOrientation = [0 0 0 0 0];

% enter manufacturer of SAT.  1 for Campbell, 0 for RMYoung.  RMYoung v = Campbell u!
info.sonicManufact = [1 1 1 1 1];  % ES5

% enter orientation of tower relative to sonic head
info.tower = 210; % ES5

% tower elevation
info.siteElevation = 1432; % (m)

% enter expected table names.  Missing tables will be filled with NaNs to create consistency 
% when multiple output files are concatnated with getData.m
info.tableNames = {'ES5_1HZ','ES5_20HZ'};

% enter table scan frequencies corresponding to tableNames
info.tableScanFrequency = [1, 20];  %[Hz]

% enter number of columns in each .csv table.  Note that the number of columns in the output structure will 
% be 3 less than the number in the .csv file.  This is because the 4 column date vector is replaced with a Matlab's 
% single-column serial time.  Also, note that View Pro frequently cuts of column 1 (the year!) of the .csv file. 
info.tableNumberOfColumns = [16 , 36];
