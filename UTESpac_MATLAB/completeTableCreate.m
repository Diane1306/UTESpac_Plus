function filledTable = completeTableCreate(beginSerialDay,endSerialDay,sampleFrequency,expectedTableColumns)

% find time step in days
dt = 1/sampleFrequency/24/60/60;

% create time stamps for filled table
filledTimeStamps = beginSerialDay + dt : dt : endSerialDay+1;

% create time stamps for filled table
% filledTimeStamps = beginSerialDay + 0 : dt : endSerialDay; % modified by Diane to use FM data start from 00:00:00 instead of 00:00:dt

% initialize filled table
filledTable = nan(numel(filledTimeStamps),expectedTableColumns);

% place campbell datevectors in columns 1:4
filledTable(:,1:4) = serialDate2CampbellDate(filledTimeStamps);
end