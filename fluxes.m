function [output, raw] = fluxes(data, rotatedSonicData, PFSonicData, info, output, sensorInfo,tableNames)
% fluxes computes all turbulent statistics.  Spike and NaN flags are used from sonic and finewires to nan-out flagged
% statistics.
fprintf('\nComputing Fluxes\n')
% ensure that sonics exist
if ~isfield(sensorInfo,'u')
    raw = [];
    return
end
try
    %---------------- FIND REFERENCE VALUES
    % constants
    Rd = 287.058;  % [J/K/kg] Gas constant for air
    Rv = 461.495;  % [J/K/kg] Gas constant for water vapor
    
    % find sonic time stamps
    t = data{1,sensorInfo.u(1,1)}(:,1); % time stamps for sonics
    
    % find reference pressure
    zRef = min(sensorInfo.u(:,3));  % zRef is lowest sonic level
    altitude = info.siteElevation;
    Pref = 101325*(1-2.25577*10^-5*(altitude+zRef))^5.25588/1000; % find pRef from elevation: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
    
    angle = info.angle; % slope angle in degrees, added by Diane

    % find pressure, use calculated Pref if barometer does not exist, is all nans or deviates from Pref by more than 5%
    if isfield(sensorInfo,'P')
        zRef = sensorInfo.P(1,3);   % Pressure Height
        Ptble = sensorInfo.P(1,1);  % Pressure table
        Pcol = sensorInfo.P(1,2);   % Pressure Column
        PkPa = data{1,Ptble}(:,Pcol);  % Pressure in kPa or mbar
        Ptime = data{1,Ptble}(:,1); % pressure time stamps
        
        % ensure that Pressure is in kPa
        if median(PkPa, 'omitmissing')>200
            PkPa = PkPa./10;
        end
        
        % store PkPa in raw tables
        if info.saveRawConditionedData
            raw.P = [Ptime, PkPa]; % u
        end
        
        % average P down to avgPer
        PkPaAvg = simpleAvg([PkPa, Ptime],info.avgPer,0);
        clear P Ptime
        if sum(isnan(PkPaAvg))/numel(PkPaAvg) < 1
            fprintf('Barometer found. Median Pressure = %0.03g kPa\n',median(PkPa, 'omitmissing'))
        end
    end
    
    if ~isfield(sensorInfo,'P') || sum(isnan(PkPaAvg))/numel(PkPaAvg) == 1 || abs((median(PkPaAvg, 'omitmissing') - Pref)/Pref) > 0.05
        zRef = min(sensorInfo.u(:,3));  % zRef is lowest sonic level
        altitude = info.siteElevation;
        Pref = 101325*(1-2.25577*10^-5*(altitude+zRef))^5.25588/1000; % find pRef from elevation: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
        PkPaAvg = Pref*ones(size(output.rotatedSonic(:,1))); % create pressure array
        fprintf('No barometer found or barometer values all NaNs or median barometer values deviate from reference pressure by more than 5%%.\nReference Pressure of %0.03g kPa calculated from site elevation (defined in siteInfo.m) will be used.\n',Pref)
    end
    
    % find reference temperature nearest to zRef
    if isfield(sensorInfo,'T')  % look for slow repsonse temperature first (HMP45/155)
        % find temperature closest to zRef
        [~, Tindex] = min(abs(sensorInfo.T(:,3)-zRef));  % find T closest to zRef
        T_tble = sensorInfo.T(Tindex,1);  % Find sonic corresponding to pressure level
        T_col = sensorInfo.T(Tindex,2);
        Tref_K = data{1,T_tble}(:,T_col); % reference sonic temperature in C or K
        Tref_time = data{1,T_tble}(:,1);
        Tref_Kavg = simpleAvg([Tref_K, Tref_time],info.avgPer,0); % average down to avgPer
        
        
        % ensure Tref is in K
        if median(Tref_Kavg, 'omitmissing')<200
            Tref_Kavg = Tref_Kavg + 273.15;  % put sonTref_K in K
        end
        
        if sum(isnan(Tref_Kavg))/numel(Tref_Kavg) < 1
            fprintf('Slow Repsonse Temperature found. Median Reference Temperature = %0.03g C\n',median(Tref_Kavg, 'omitmissing')-273.15)
        end
    end
    
    if ~isfield(sensorInfo,'T') || sum(isnan(Tref_Kavg))/numel(Tref_Kavg) == 1 % if slow response T unavailable, uses closest sonic to zRef
        % find sonic reference temperature closest to zRef
        [~, Tindex] = min(abs(sensorInfo.u(:,3)-zRef));  % find sonic closest to zRef
        T_tble = sensorInfo.Tson(Tindex,1);
        T_col = sensorInfo.Tson(Tindex,2);
        Tref_K = data{1,T_tble}(:,T_col); % reference sonic temperature in C or K
        
        Tref_Kavg = simpleAvg([Tref_K, t],info.avgPer,0); % average down to avgPer
        
        
        % ensure Tref is in K
        if median(Tref_Kavg, 'omitmissing')<200
            Tref_Kavg = Tref_Kavg + 273.15;  % put sonTref_K in K
        end
        fprintf('No Slow Repsonse Temperature found, or slow response temperature is all NaNs! Median Reference Temperature from Sonic = %0.03g C\n',median(Tref_Kavg, 'omitmissing')-273.15)
    end
    
    % find reference specific humidity, this will be updated at each level with the nearest HMP if available
    if isfield(sensorInfo,'RH')  % try slow response RH sensor first
        
        % find closes RH measurement to zRef
        [~, RHindex] = min(abs(sensorInfo.RH(:,3)-zRef));
        RHtble = sensorInfo.RH(RHindex,1);
        RHcol = sensorInfo.RH(RHindex,2);
        RH = data{1,RHtble}(:,RHcol);  % RH closest to zRef
        RHtime = data{1,RHtble}(:,1);
        
        % interpolate RH to same time stamps as sonic_time
        RHavg = simpleAvg([RH, RHtime],info.avgPer,1);
        RHavg_t = RHavg(:,2);
        RHavg(:,2) = [];
        clear RH
        
        % find average specific humidity
        qRefavg = RHtoSpecHum(RHavg,PkPaAvg,Tref_Kavg); % [kg/kg]
        
        % ensure that qRefavg is not all NaNs
        if sum(isnan(qRefavg))<numel(qRefavg)
            % check output here: http://www.rotronic.com/humidity_measurement-feuchtemessung-mesure_de_l_humidite/humidity-calculator-feuchterechner-mr
            % interpolate to sonic time stamps.  Pad qRefavg at beginning to allow interpolation at all time stamps
            x = RHavg_t(~isnan(qRefavg));
            x = [floor(x(1)); x];
            y = qRefavg(~isnan(qRefavg));
            y = [y(1); y];
            qRefFast = interp1(x,y,t);  % interpolate qRef to Sonic Frequency
            fprintf('Slow Repsonse RH found! Median Reference Specific Humidity = %0.03g g/kg\n',1000*median(qRefavg, "omitmissing"))
        end
    else
        qRefavg = nan;
    end
    
    % find qref from IRGA or KH2O
    if ~isfield(sensorInfo,'RH') && (isfield(sensorInfo,'irgaH2O') || isfield(sensorInfo,'KH2O')) 
        
        if isfield(sensorInfo,'irgaH2O')
            rhoH2O_tble = sensorInfo.irgaH2O(1,1);
            rhoH2O_col = sensorInfo.irgaH2O(1,2);
            rhovIRGA = data{1,rhoH2O_tble}(:,rhoH2O_col);
            rhov_t = data{1,rhoH2O_tble}(:,1);
            rhovIRGAavg = simpleAvg([rhovIRGA,rhov_t],info.avgPer);
            disp('No slow-response RH found. qRef calculated from EC150')
        elseif isfield(sensorInfo,'KH2O')
            rhoH2O_tble = sensorInfo.KH2O(1,1);
            rhoH2O_col = sensorInfo.KH2O(1,2);
            rhovIRGA = data{1,rhoH2O_tble}(:,rhoH2O_col);
            rhov_t = data{1,rhoH2O_tble}(:,1);
            rhovIRGAavg = simpleAvg([rhovIRGA,rhov_t],info.avgPer);
            disp('No slow-response RH found. qRef calculated from KH2O')
        elseif isfield(sensorInfo,'LiH2O')
            rhoH2O_tble = sensorInfo.LiH2O(1,1);
            rhoH2O_col = sensorInfo.LiH2O(1,2);
            rhovIRGA = data{1,rhoH2O_tble}(:,rhoH2O_col)*18/1000;  % Multiply by 18/1000 to go from mmol/mol to g/m^3
            rhov_t = data{1,rhoH2O_tble}(:,1);
            rhovIRGAavg = simpleAvg([rhovIRGA,rhov_t],info.avgPer);
            disp('No slow-response RH found. qRef calculated from EC150') 
        end
        qRefavg = (rhovIRGAavg(:,1)./1000)./(PkPaAvg*1000./(Rd*Tref_Kavg)); % kg/kg  rhoH2O/rhoAir
        x = rhovIRGAavg(~isnan(qRefavg),2);
        x = [floor(x(1)); x];
        y = qRefavg(~isnan(qRefavg));
        y = [y(1); y];
        qRefFast = interp1(x,y,t);  % interpolate qRef to Sonic Frequency
        fprintf('Median Reference Specific Humidity = %0.03g g/kg\n',1000*median(qRefavg, 'omitmissing'))
    elseif ~isfield(sensorInfo,'RH') ||  sum(isnan(qRefavg))/numel(qRefavg)==1 % if no humidity measurements exist, use info.qRef
        qRefavg = info.qRef./1000*ones(size(Tref_Kavg));
        qRefFast = info.qRef./1000*ones(size(t));
        fprintf('No slow-response nor fast response RH found (or humidity measurement is all NaNs). qRef = %0.02g g/kg defined in INORMATION block is being used.',info.qRef)
    end
    % store specific humidity in output
    
    % find moist-air density, dry-air density and vapor density
    PvAvg = qRefavg.*PkPaAvg./0.622;  % partial pressure of water vapor (kPa)
    PdAvg = PkPaAvg-PvAvg; % partial pressure of dry air
    rhodAvg = 1000*PdAvg./(Rd*Tref_Kavg); % density of dry air (kg/m^3)
    rhovAvg = 1000*PvAvg./(Rv*Tref_Kavg); % density of water vapor (kg/m^3)
    rhoAvg = rhodAvg+rhovAvg; % density of moist air (kg/m^3)
    fprintf('Calculated Moist Air Density = %0.03g kg/m^3, Dry Air Density = %0.03g kg/m^3\n',median(rhoAvg, 'omitmissing'),median(rhodAvg, 'omitmissing'))
    
    
    %------------- ITERATE THROUGH SONICS
    % find number of sonics
    numSonics = size(sensorInfo.u,1);
    
    for ii = 1:numSonics
        if ii==1
            
            % find total number of averaging periods
            N = round((t(end)-t(1))/(info.avgPer/(24*60)));
            
            % find breakpoints for detrending.  size(bp,1) = size(N,1) + 1
            bp = round(linspace(0,numel(t),N+1));
            
            % initialize flux matrices
            H = nan(N,1);  % kinematic sensible heat flux
            Hlat = nan(N,1); % kinematic sensible, lateral heat flux
            tau = nan(N,1); % momentum flux
            tke = nan(N,1); % turbulent kinetic energy
            LHflux = nan(N,1); % latent heat flux
            CO2flux = nan(N,1); % CO2 flux
            derivedT = nan(N,1);  % matrix for derived temperatures
            L = nan(N,1); % Obukhov Length
            sigma = nan(N,1); % standard deviations (u, v, w)

            % initialize correlation coefficient matrices by Diane
            % 2025/08/05
            R = nan(N,1);

            % initialize flux transport efficiencies by Diane
            eta = nan(N,1);
            delta_flux_ctrb = nan(N,1);
            delta_time_ctrb = nan(N,1);

            % initialize turubulent transport and dissipation terms for TKE
            % budget
            turbtr = nan(N,1);
            epsilon = nan(N,1);

            % initialize skewness by Diane
            skew = nan(N,1);

            % initialize Heat flux in SNSP system
            H_SNSP = nan(N,1);

            % initialize Lateral fluxes for all scalars
            Flux_lat = nan(N,1);
            
            % initialize raw variable matrices
            if info.saveRawConditionedData
%                 raw.u = nan(bp(end),numSonics); % u
%                 raw.v = nan(bp(end),numSonics); % v
%                 raw.w = nan(bp(end),numSonics); % w
%                 raw.WD = nan(bp(end),numSonics);  % raw wind direction
                raw.uPF = nan(bp(end),numSonics);  % u planar fit and yaw corrected
                raw.vPF = nan(bp(end),numSonics);  % v planar fit and yaw corrected
                raw.wPF = nan(bp(end),numSonics);  % w planar fit and yaw corrected
                raw.u_tilt = nan(bp(end),numSonics);  % u planar fit cross slope
                raw.v_tilt = nan(bp(end),numSonics);  % v planar fit along slope
                raw.w_tilt = nan(bp(end),numSonics);  % w planar fit slope normal
                raw.sonTs = nan(bp(end),numSonics); % Temperature from Sonic
                raw.Theta_v_son = nan(bp(end),numSonics); % virtual potential Temperature from Sonic
%                 raw.uPF_Prime = nan(bp(end),numSonics);  % u' planar fit and yaw corrected
%                 raw.vPF_Prime = nan(bp(end),numSonics); % v' planar fit and yaw corrected
%                 raw.wPF_Prime = nan(bp(end),numSonics); % w' planar fit and yaw corrected
%                 raw.u_tilt_Prime = nan(bp(end),numSonics);  % u' planar fit 
%                 raw.v_tilt_Prime = nan(bp(end),numSonics); % v' planar fit 
%                 raw.w_tilt_Prime = nan(bp(end),numSonics); % w' planar fit 
%                 raw.sonTsPrime = nan(bp(end),numSonics); % Ts' from sonic
%                 raw.Theta_v_sonPrime = nan(bp(end),numSonics);
                raw.t = t; % serial time stamp
                raw.z = nan(1,numSonics); % sonic heights
                
                % finewires
                if isfield(sensorInfo,'fw')
                    raw.fwTh = nan(bp(end),numSonics); % theta from FW
                    raw.fwT = nan(bp(end),numSonics); % T from FW
                    raw.fwTPrime = nan(bp(end),numSonics); % T' from FW
                    raw.fwThPrime = nan(bp(end),numSonics); % theta' from finewire
                end
                
                
                % water vapor
                if isfield(sensorInfo,'irgaH2O')
                    numH2OSensors = size(sensorInfo.irgaH2O,1);
                    raw.rhov = nan(bp(end),numH2OSensors); % H2O
                    raw.rhovPrime = nan(bp(end),numH2OSensors); % H2O'
                elseif isfield(sensorInfo,'LiH2O')
                    numH2OSensors = size(sensorInfo.LiH2O,1);
                    raw.rhov = nan(bp(end),numH2OSensors); % H2O
                    raw.rhovPrime = nan(bp(end),numH2OSensors); % H2O'
                elseif isfield(sensorInfo,'KH2O')
                    numH2OSensors = size(sensorInfo.KH2O,1);
                    raw.rhov = nan(bp(end),numH2OSensors); % H2O
                    raw.rhovPrime = nan(bp(end),numH2OSensors); % H2O'
                end
                
                % CO2
                if isfield(sensorInfo,'irgaCO2')
                    numCO2Sensors = size(sensorInfo.irgaCO2,1);
                    raw.rhoCO2 = nan(bp(end),numCO2Sensors); % CO2
                    raw.rhoCO2Prime = nan(bp(end),numCO2Sensors); % CO2'
                    raw.rhoCO2exteralPrime = nan(bp(end),numCO2Sensors); % external CO2'
                elseif isfield(sensorInfo,'LiCO2')
                    numCO2Sensors = size(sensorInfo.LiCO2,1);
                    raw.rhoCO2 = nan(bp(end),numCO2Sensors); % CO2
                    raw.rhoCO2Prime = nan(bp(end),numCO2Sensors); % CO2'
                    raw.rhoCO2exteralPrime = nan(bp(end),numCO2Sensors); % external CO2'
                end
                
            else
                raw = [];
            end
            
            % initialize headers
            Hheader = cell(1);
            HlatHeader = cell(1);
            tauHeader = cell(1);
            tkeHeader = cell(1);
            LHfluxHeader = cell(1);
            CO2fluxHeader = cell(1);
            derivedTheader = cell(1);
            Lheader = cell(1);
            sigmaHeader = cell(1);
            RHeader = cell(1);
            etaHeader = cell(1);
            delta_flux_ctrbHeader = cell(1);
            delta_time_ctrbHeader = cell(1);
            turbtrHeader = cell(1);
            epsilonHeader = cell(1);
            skewHeader = cell(1);
            H_SNSPHeader = cell(1);
            Flux_latHeader = cell(1);
        end
        try
            
            % find sonic information
            tble = sensorInfo.u(ii,1);
            sonHeight = sensorInfo.u(ii,3);
            uCol = sensorInfo.u(sensorInfo.u(:,3)==sonHeight,2);
            vCol = sensorInfo.v(sensorInfo.v(:,3)==sonHeight,2);
            wCol = sensorInfo.w(sensorInfo.v(:,3)==sonHeight,2);
            TsCol = sensorInfo.Tson(sensorInfo.Tson(:,3)==sonHeight,2);
            
            % find sonic flag information (averaged values!)
            nanFlagTableName = [tableNames{tble},'NanFlag'];
            spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
            uNanFlag = output.(nanFlagTableName)(:,uCol);
            uSpikeFlag = output.(spikeFlagTableName)(:,uCol);
            vNanFlag = output.(nanFlagTableName)(:,vCol);
            vSpikeFlag = output.(spikeFlagTableName)(:,vCol);
            wNanFlag = output.(nanFlagTableName)(:,wCol);
            wSpikeFlag = output.(spikeFlagTableName)(:,wCol);
            TsNanFlag = output.(nanFlagTableName)(:,TsCol);
            TsSpikeFlag = output.(nanFlagTableName)(:,TsCol);
            
            % check for diagnostic
            if isfield(sensorInfo,'sonDiagnostic')
                sonicDiagnosticCol = sensorInfo.sonDiagnostic(sensorInfo.sonDiagnostic(:,3)==sonHeight,2);
                if isempty(sonicDiagnosticCol)
                    sonicDiagnosticFlag = zeros(size(uNanFlag));
                else
                    sonicDiagnosticFlag = output.(tableNames{tble})(:,sonicDiagnosticCol);
                    sonicDiagnosticFlag(sonicDiagnosticFlag < info.diagnosticTest.meanSonicDiagnosticLimit) = 0;
                end
            else
                sonicDiagnosticFlag = zeros(size(uNanFlag));
            end
            sonicDiagnosticFlag(isnan(sonicDiagnosticFlag)) = 0;
            
            % find unrotated sonic values
            u = data{1,tble}(:,uCol);
            v = data{1,tble}(:,vCol);
            w = data{1,tble}(:,wCol);
            unrotatedSonFlag = logical(wNanFlag+wSpikeFlag+sonicDiagnosticFlag); % total sonic flag for unrotated calculations
            
            % find rotated sonic columns
            uCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:u'));
            vCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:v'));
            wCol = strcmp(output.rotatedSonicHeader,strcat(num2str(sonHeight),'m:w'));
            
            % load rotated sonic values
            uPF = rotatedSonicData(:,uCol);
            vPF = rotatedSonicData(:,vCol);
            wPF = rotatedSonicData(:,wCol);
            rotatedSonFlag = logical(uNanFlag+uSpikeFlag+vNanFlag+vSpikeFlag+wNanFlag+wSpikeFlag+sonicDiagnosticFlag); % total sonic flag for rotated calculations
            
            %load tilted sonic values in sonic coordinate system
            u_tilt = PFSonicData(:,uCol);
            v_tilt = PFSonicData(:,vCol);
            w_tilt = PFSonicData(:,wCol);
            % TiltedSonFlag = rotatedSonFlag;

            % find wind direction
            DirCol = strcmp(output.spdAndDirHeader,strcat(num2str(sonHeight),'m direction'));
            direction = output.spdAndDir(:, DirCol);

            % find sonic temperature
            Tson = data{1,tble}(:,TsCol);
            if median(Tson, 'omitmissing') > 250  % ensure temperature is in C
                Tson = Tson - 273.15;
            end
            TsonFlag = logical(TsNanFlag+TsSpikeFlag);
            
            % find fw
            if isfield(sensorInfo,'fw')
                fwCol = sensorInfo.fw(sensorInfo.fw(:,3)==sonHeight,2);
                fw = data{1,tble}(:,fwCol);
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                fwNanFlag = output.(nanFlagTableName)(:,fwCol);
                fwSpikeFlag = output.(spikeFlagTableName)(:,fwCol);
                fwFlag = logical(fwNanFlag+fwSpikeFlag);
            else
                fw = [];
                fwFlag = [];
            end
            
            % find specific humidity at level if HMP exists
            if isfield(sensorInfo,'RH') && isfield(sensorInfo,'T') && any(sensorInfo.RH(:,3) == sonHeight)
                
                % find RH at sensor level
                RHtble = sensorInfo.RH(sensorInfo.RH(:,3)==sonHeight,1);
                RHcol = sensorInfo.RH(sensorInfo.RH(:,3)==sonHeight,2);
                RHavg = simpleAvg([data{1,RHtble}(:,RHcol), data{1,RHtble}(:,1)],info.avgPer,0);
                
                % find T at sensor level
                Ttble = sensorInfo.T(sensorInfo.T(:,3)==sonHeight,1);
                Tcol = sensorInfo.T(sensorInfo.T(:,3)==sonHeight,2);
                Tavg = simpleAvg([data{1,Ttble}(:,Tcol), data{1,Ttble}(:,1)],info.avgPer,0);
                
                % put Tavg in K
                if median(Tavg, 'omitmissing') < 250; Tavg = Tavg + 273.15; end
                
                % intialize specific humidity output, populate col 1 with time stamps
                if ~isfield(output,'specificHum')
                    output.specificHum = t(bp(2:end));
                    output.specificHumHeader = cell(1);
                    output.specificHumHeader{1} = 'time';
                end
                
                % find qRefLocal and qRefFastLocal
                qRefavgLocal = RHtoSpecHum(RHavg,PkPaAvg,Tavg);
                x = RHavg_t(~isnan(qRefavgLocal));
                if isempty(x)  % if all NaNs, use qRef (non-local)
                    qRefFastLocal = qRefFast;
                    
                    % store nan'd data at height
                    output.specificHum(:,end+1) = qRefavgLocal*nan;
                    output.specificHumHeader{end+1} = sprintf('%g m: q(g/g)',sonHeight);
                else
                    x = [floor(x(1)); x];
                    y = qRefavgLocal(~isnan(qRefavgLocal));
                    y = [y(1); y];
                    qRefFastLocal = interp1(x,y,t);  % interpolate qRef to Sonic Frequency
                    
                    % store data
                    output.specificHum(:,end+1) = qRefavgLocal;
                    output.specificHumHeader{end+1} = sprintf('%g m: q(g/g)',sonHeight);
                end
                
                % use qRefFastLocal if not all NaNs to find virt temp
                if sum(isnan(qRefFastLocal))/numel(qRefFastLocal) == 1
                    qRefFastLocal = qRefFast;
                end
            else
                qRefFastLocal = qRefFast;
            end
            
            % find fw pot temp
            Gamma = 0.0098; % Dry Lapse Rate. K/m
            if ~isempty(fw)
                thetaFw = fw + Gamma*(sonHeight - zRef);
                temp = simpleAvg([thetaFw t],info.avgPer);
                derivedT(:,end+1) = temp(:,1);
                derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_fw');
            else
                thetaFw = [];
            end
            
            % find sonic pot temp
            % thetaSon = Tson + Gamma*(sonHeight - zRef);
            thetaSon = ((Tson + 273.15)./(1+0.51*qRefFastLocal) + Gamma*(sonHeight - zRef)) .*(1+0.61*qRefFastLocal) - 273.15; % Modified by Diane to get sonic virtual potential temperature
            temp = simpleAvg([thetaSon t],info.avgPer);
            derivedT(:,end+1) = temp(:,1);
            % derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_s_son');
            derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_v_son'); % by Diane
            
            % find fw virt, pot temp with qRef
            if ~isempty(thetaFw)
                VthetaFw = thetaFw.*(1+0.61*qRefFastLocal); % stull pg 7
                temp = simpleAvg([VthetaFw t],info.avgPer);
                derivedT(:,end+1) = temp(:,1);
                derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_v_fw');
            else
                VthetaFw = [];
            end
            
            % find pot temp from sonic
            % thetaSonAir = thetaSon./(1+0.51*qRefFastLocal); % stull pg 7
            thetaSonAir = (Tson + 273.15)./(1+0.51*qRefFastLocal) - 273.15; % IRGAOSN manual p43
            temp = simpleAvg([thetaSonAir t],info.avgPer);
            derivedT(:,end+1) = temp(:,1);
            % derivedTheader{end+1} = strcat(num2str(sonHeight),' m: theta_son');
            derivedTheader{end+1} = strcat(num2str(sonHeight),' m: T_son_air'); % by Diane
            
            % find H2O and CO2 columns if they exist
            if isfield(sensorInfo,'irgaH2O') && ~isempty(sensorInfo.irgaH2O(sensorInfo.irgaH2O(:,3)==sonHeight,2)) % EC150
                irgaH2Ocol = sensorInfo.irgaH2O(sensorInfo.irgaH2O(:,3)==sonHeight,2);
                irgaGasDiagCol = sensorInfo.irgaGasDiag(sensorInfo.irgaGasDiag(:,3)==sonHeight,2);
                irgaH2OSigCol = sensorInfo.irgaH2OsigStrength(sensorInfo.irgaH2OsigStrength(:,3)==sonHeight,2);
                irgaH2O = data{1,tble}(:,irgaH2Ocol);
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                H2ONanFlag = output.(nanFlagTableName)(:,irgaH2Ocol);
                H2OSpikeFlag = output.(spikeFlagTableName)(:,irgaH2Ocol);
                H2OsigFlag = output.(tableNames{tble})(:,irgaH2OSigCol);
                H2OsigFlag(H2OsigFlag > info.diagnosticTest.H2OminSignal) = 0;
                H2OsigFlag(isnan(H2OsigFlag)) = 0;
                gasDiagFlag = output.(tableNames{tble})(:,irgaGasDiagCol);
                gasDiagFlag(gasDiagFlag < info.diagnosticTest.meanGasDiagnosticLimit) = 0;
                gasDiagFlag(isnan(gasDiagFlag)) = 0;
                H2OFlag = logical(H2ONanFlag+H2OSpikeFlag+H2OsigFlag+gasDiagFlag);
            elseif isfield(sensorInfo,'LiH2O') && ~isempty(sensorInfo.LiH2O(sensorInfo.LiH2O(:,3)==sonHeight,2)) % Li7500
                irgaH2Ocol = sensorInfo.LiH2O(sensorInfo.LiH2O(:,3)==sonHeight,2);
                irgaH2O = data{1,tble}(:,irgaH2Ocol)*0.018; % convert from mmol/m^3 to g/m^3
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                H2ONanFlag = output.(nanFlagTableName)(:,irgaH2Ocol);
                H2OSpikeFlag = output.(spikeFlagTableName)(:,irgaH2Ocol);
                % check for LiCor diagnostic flag
                if isfield(sensorInfo,'LiGasDiag')
                    irgaGasDiagCol = sensorInfo.LiGasDiag(sensorInfo.LiGasDiag(:,3)==sonHeight,2);
                    gasDiagFlag = output.(tableNames{tble})(:,irgaGasDiagCol);
                    gasDiagFlag(gasDiagFlag > info.diagnosticTest.meanLiGasDiagnosticLimit) = 0;
                    gasDiagFlag(isnan(gasDiagFlag)) = 0;
                else
                    gasDiagFlag = false(size(H2ONanFlag));                    
                end               
                H2OFlag = logical(H2ONanFlag+H2OSpikeFlag+gasDiagFlag);
            else
                irgaH2O = [];
                H2OFlag = [];
            end
            
            if isfield(sensorInfo,'irgaCO2') && ~isempty(sensorInfo.irgaCO2(sensorInfo.irgaCO2(:,3)==sonHeight,2)) % EC150
                irgaCO2col = sensorInfo.irgaCO2(sensorInfo.irgaCO2(:,3)==sonHeight,2);
                irgaCO2SigCol = sensorInfo.irgaCO2sigStrength(sensorInfo.irgaCO2sigStrength(:,3)==sonHeight,2);
                irgaCO2 = data{1,tble}(:,irgaCO2col);  
                CO2sigFlag = output.(tableNames{tble})(:,irgaCO2SigCol);
                CO2sigFlag(CO2sigFlag > info.diagnosticTest.CO2minSignal) = 0;
                CO2sigFlag(isnan(CO2sigFlag)) = 0;
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                CO2NanFlag = output.(nanFlagTableName)(:,irgaCO2col);
                CO2SpikeFlag = output.(spikeFlagTableName)(:,irgaCO2col);
                CO2Flag = logical(CO2NanFlag+CO2SpikeFlag+gasDiagFlag+CO2sigFlag);
            elseif isfield(sensorInfo,'LiCO2') && ~isempty(sensorInfo.LiCO2(sensorInfo.LiCO2(:,3)==sonHeight,2))%Li7500
                irgaCO2col = sensorInfo.LiCO2(sensorInfo.LiCO2(:,3)==sonHeight,2);
                irgaCO2 = data{1,tble}(:,irgaCO2col)*44; % Mult by 44 to convert from mmol/m^3 to mg/m^3
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                CO2NanFlag = output.(nanFlagTableName)(:,irgaCO2col);
                CO2SpikeFlag = output.(spikeFlagTableName)(:,irgaCO2col);
                CO2Flag = logical(CO2NanFlag+CO2SpikeFlag+gasDiagFlag);
            else
                irgaCO2 = [];
                CO2Flag = [];
            end
            if isfield(sensorInfo,'KH2O')
                KH2Ocol = sensorInfo.KH2O(sensorInfo.KH2O(:,3)==sonHeight,2);
                KH2O = data{1,tble}(:,KH2Ocol);
                nanFlagTableName = [tableNames{tble},'NanFlag'];
                spikeFlagTableName = [tableNames{tble},'SpikeFlag'];
                H2ONanFlag = output.(nanFlagTableName)(:,KH2Ocol);
                H2OSpikeFlag = output.(spikeFlagTableName)(:,KH2Ocol);
                H2OFlag = logical(H2ONanFlag+H2OSpikeFlag);
            elseif isfield(sensorInfo,'irgaH2O') || isfield(sensorInfo,'LiH2O') % don't delete H2O flag if irga exists!
                KH2O = [];
            else
                KH2O = [];
                H2OFlag = [];
            end
            
            %--------------------- ITERATE THROUGH ALL TIME STEPS
            for jj = 1:N
                if ii == 1
                    % place time stamp in column 1
                    H(jj,1) = t(bp(jj+1));
                    Hlat(jj,1) = t(bp(jj+1));
                    tau(jj,1) = t(bp(jj+1));
                    tke(jj,1) = t(bp(jj+1));
                    LHflux(jj,1) = t(bp(jj+1));
                    CO2flux(jj,1) = t(bp(jj+1));
                    derivedT(jj,1) = t(bp(jj+1));
                    L(jj,1) = t(bp(jj+1));
                    sigma(jj,1) = t(bp(jj+1));
                    R(jj,1) = t(bp(jj+1));
                    eta(jj,1) = t(bp(jj+1));
                    delta_flux_ctrb(jj,1) = t(bp(jj+1));
                    delta_time_ctrb(jj,1) = t(bp(jj+1));
                    turbtr(jj,1) = t(bp(jj+1));
                    epsilon(jj,1) = t(bp(jj+1));
                    skew(jj,1) = t(bp(jj+1));
                    H_SNSP(jj,1) = t(bp(jj+1));
                    Flux_lat(jj,1) = t(bp(jj+1));
                    if jj == 1
                        Hheader{1} = 'time';
                        HlatHeader{1} = 'time';
                        tauHeader{1} = 'time';
                        tkeHeader{1} = 'time';
                        derivedTheader{1} = 'time';
                        LHfluxHeader{1} = 'time';
                        CO2fluxHeader{1} = 'time';
                        Lheader{1} = 'time';
                        sigmaHeader{1} = 'time';
                        RHeader{1} = 'time';
                        etaHeader{1} = 'time';
                        delta_flux_ctrbHeader{1} = 'time';
                        delta_time_ctrbHeader{1} = 'time';
                        turbtrHeader{1} = 'time';
                        epsilonHeader{1} = 'time';
                        skewHeader{1} = 'time';
                        H_SNSPHeader{1} = 'time';
                        Flux_latHeader{1} = 'time';
                    end
                end
                
                % place rho in column 2 of H
                Hheader{2} = 'rho';
                H(jj,2) = rhoAvg(jj);
                
                % place Cp in column 3 of H
                Hheader{3} = 'cp';
                H(jj,3) = 1004.67*(1+0.84*qRefavg(jj)); % pg. 640 of Stull
                
                % find unrotated sonic pertubations
                uP = nandetrend(u(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                vP = nandetrend(v(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                wP = nandetrend(w(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                TsonP = nandetrend(Tson(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                % Theta_v_sonP = nandetrend(thetaSonAir(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                Theta_v_sonP = nandetrend(thetaSon(bp(jj)+1:bp(jj+1)),info.detrendingFormat); % Modified by Diane to get sonic virtual potential temperature
                T_air_sonP = nandetrend(thetaSonAir(bp(jj)+1:bp(jj+1)),info.detrendingFormat); % Added by Diane to get air temperature from sonic measurements
                T_air_sonM = mean(thetaSonAir(bp(jj)+1:bp(jj+1)),'omitmissing'); % Added by Diane to get air temperature from sonic measurements
                
                % find rotated sonic pertubations
                uPF_P = nandetrend(uPF(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                vPF_P = nandetrend(vPF(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                wPF_P = nandetrend(wPF(bp(jj)+1:bp(jj+1)),info.detrendingFormat);

                % find tilted sonic perturbations
                u_tilt_P = nandetrend(v_tilt(bp(jj)+1:bp(jj+1)),info.detrendingFormat); % v_tilt positive as downslope
                v_tilt_P = nandetrend(u_tilt(bp(jj)+1:bp(jj+1)),info.detrendingFormat); % u_tilt cross-slope
                w_tilt_P = nandetrend(w_tilt(bp(jj)+1:bp(jj+1)),info.detrendingFormat); % w_tilt slope-normal
                
                % store rotated data
                if info.saveRawConditionedData
%                     raw.u(bp(jj)+1:bp(jj+1),ii) = u(bp(jj)+1:bp(jj+1));
%                     raw.v(bp(jj)+1:bp(jj+1),ii) = v(bp(jj)+1:bp(jj+1));
%                     raw.w(bp(jj)+1:bp(jj+1),ii) = w(bp(jj)+1:bp(jj+1));
                    raw.uPF(bp(jj)+1:bp(jj+1),ii) = uPF(bp(jj)+1:bp(jj+1));
                    raw.vPF(bp(jj)+1:bp(jj+1),ii) = vPF(bp(jj)+1:bp(jj+1));
                    raw.wPF(bp(jj)+1:bp(jj+1),ii) = wPF(bp(jj)+1:bp(jj+1));
                    raw.u_tilt(bp(jj)+1:bp(jj+1),ii) = u_tilt(bp(jj)+1:bp(jj+1));
                    raw.v_tilt(bp(jj)+1:bp(jj+1),ii) = v_tilt(bp(jj)+1:bp(jj+1));
                    raw.w_tilt(bp(jj)+1:bp(jj+1),ii) = w_tilt(bp(jj)+1:bp(jj+1));
%                     raw.uPF_Prime(bp(jj)+1:bp(jj+1),ii) = uPF_P;
%                     raw.vPF_Prime(bp(jj)+1:bp(jj+1),ii) = vPF_P;
%                     raw.wPF_Prime(bp(jj)+1:bp(jj+1),ii) = wPF_P;
%                     raw.u_tilt_Prime(bp(jj)+1:bp(jj+1),ii) = u_tilt_P;
%                     raw.v_tilt_Prime(bp(jj)+1:bp(jj+1),ii) = v_tilt_P;
%                     raw.w_tilt_Prime(bp(jj)+1:bp(jj+1),ii) = w_tilt_P;
                    raw.z(ii) = sonHeight;
                end
                
                % find standard deviations of wind vector, sonic and finewire temperature
                if ~isempty(fw) % check for fw first
                    numSigmaVariables = 9;
                    
                    sigma(jj,10+(ii-1)*numSigmaVariables) = std(fw(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{10+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_TFW');
                    if fwFlag(jj); sigma(jj,10+(ii-1)*numSigmaVariables) = nan; end;
                    
                else
                    % numSigmaVariables = 8;
                    numSigmaVariables = 12; % Diane added 4 columns for theta_v, q, CO2, CO2_WPL
                end
                
                sigma(jj,2+(ii-1)*numSigmaVariables) = std(u(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{2+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_u');
                if rotatedSonFlag(jj); sigma(jj,2+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,3+(ii-1)*numSigmaVariables) = std(v(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{3+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_v');
                if rotatedSonFlag(jj); sigma(jj,3+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,4+(ii-1)*numSigmaVariables) = std(w(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{4+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_w');
                if unrotatedSonFlag(jj); sigma(jj,4+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,5+(ii-1)*numSigmaVariables) = std(uPF(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{5+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_uPF');
                if rotatedSonFlag(jj); sigma(jj,5+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,6+(ii-1)*numSigmaVariables) = std(vPF(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{6+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_vPF');
                if rotatedSonFlag(jj); sigma(jj,6+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,7+(ii-1)*numSigmaVariables) = std(wPF(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{7+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_wPF');
                if rotatedSonFlag(jj); sigma(jj,7+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,8+(ii-1)*numSigmaVariables) = std(Tson(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{8+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_Tson');
                if TsonFlag(jj); sigma(jj,8+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,9+(ii-1)*numSigmaVariables) = mean(wPF_P.*TsonP.^2, 'omitmissing'); sigmaHeader{9+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :wPFP_TsonP_TsonP');
                if TsonFlag(jj) || rotatedSonFlag(jj); sigma(jj,9+(ii-1)*numSigmaVariables) = nan; end;
                
                sigma(jj,10+(ii-1)*numSigmaVariables) = std(thetaSon(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{10+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_Theta_v');
                if TsonFlag(jj) || rotatedSonFlag(jj); sigma(jj,10+(ii-1)*numSigmaVariables) = nan; end;
                
                % find correlation coefficients between fluxes
                numRVariables = 10; % 5 for fluxes, 5 between quantities (including CO2_WPL)
                fuw = uPF_P.*wPF_P; fwt = Theta_v_sonP.*wPF_P; indx = ~(isnan(fuw) | isnan(fwt));
                R_mat_temp = corrcoef(fuw(indx), fwt(indx));
                R(jj, 2+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{2+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_uPFwPF_wPFThetav');
                if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 2+(ii-1)*numRVariables) = nan; end;
                
                % find transport efficiencies for momentum and heat
                numetaVariables = 4;
                eta(jj, 2+(ii-1)*numetaVariables) = find_eta(wPF_P, uPF_P); etaHeader{2+(ii-1)*numetaVariables} = strcat(num2str(sonHeight),'m :eta_wPFuPF');
                if TsonFlag(jj) || rotatedSonFlag(jj); eta(jj, 2+(ii-1)*numetaVariables) = nan; end;
                
                eta(jj, 3+(ii-1)*numetaVariables) = find_eta(wPF_P, Theta_v_sonP); etaHeader{3+(ii-1)*numetaVariables} = strcat(num2str(sonHeight),'m :eta_wPFThetav');
                if TsonFlag(jj) || rotatedSonFlag(jj); eta(jj, 3+(ii-1)*numetaVariables) = nan; end;
                

                numDelta_flux_ctrbVariables = 4;
                delta_flux_ctrb(jj, 2+(ii-1)*numDelta_flux_ctrbVariables) = find_delta_flux(wPF_P, uPF_P); delta_flux_ctrbHeader{2+(ii-1)*numDelta_flux_ctrbVariables} = strcat(num2str(sonHeight),'m :S_wPFuPF');
                if TsonFlag(jj) || rotatedSonFlag(jj); delta_flux_ctrb(jj, 2+(ii-1)*numDelta_flux_ctrbVariables) = nan; end;
                
                delta_flux_ctrb(jj, 3+(ii-1)*numDelta_flux_ctrbVariables) = find_delta_flux(wPF_P, Theta_v_sonP); delta_flux_ctrbHeader{3+(ii-1)*numDelta_flux_ctrbVariables} = strcat(num2str(sonHeight),'m :S_wPFThetav');
                if TsonFlag(jj) || rotatedSonFlag(jj); delta_flux_ctrb(jj, 3+(ii-1)*numDelta_flux_ctrbVariables) = nan; end;
                

                numDelta_time_ctrbVariables = 4;
                delta_time_ctrb(jj, 2+(ii-1)*numDelta_time_ctrbVariables) = find_delta_time(wPF_P, uPF_P); delta_time_ctrbHeader{2+(ii-1)*numDelta_time_ctrbVariables} = strcat(num2str(sonHeight),'m :D_wPFuPF');
                if TsonFlag(jj) || rotatedSonFlag(jj); delta_time_ctrb(jj, 2+(ii-1)*numDelta_time_ctrbVariables) = nan; end;
                
                delta_time_ctrb(jj, 3+(ii-1)*numDelta_time_ctrbVariables) = find_delta_time(wPF_P, Theta_v_sonP); delta_time_ctrbHeader{3+(ii-1)*numDelta_time_ctrbVariables} = strcat(num2str(sonHeight),'m :D_wPFThetav');
                if TsonFlag(jj) || rotatedSonFlag(jj); delta_time_ctrb(jj, 3+(ii-1)*numDelta_time_ctrbVariables) = nan; end;
                

                % find rotated and unrated momentum flux
                numTauVariables = 14;
                tau(jj,2+(ii-1)*numTauVariables) = sqrt(mean(uP.*wP, 'omitmissing')^2+mean(vP.*wP, 'omitmissing')^2); tauHeader{2+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),'m :sqrt(u''w''^2+v''w''^2)');
                if rotatedSonFlag(jj); tau(jj,2+(ii-1)*numTauVariables) = nan; end;
                
                tau(jj,3+(ii-1)*numTauVariables) = sqrt(mean(uPF_P.*wPF_P, 'omitmissing')^2+mean(vPF_P.*wPF_P, 'omitmissing')^2); tauHeader{3+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),'m :sqrt(uPF''wPF''^2+vPF''wPF''^2)');
                if rotatedSonFlag(jj); tau(jj,3+(ii-1)*numTauVariables) = nan; end;
                
                tau(jj,4+(ii-1)*numTauVariables) = mean(uPF_P.*uPF_P, 'omitmissing'); tauHeader{4+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :uPF'uPF'");
                if rotatedSonFlag(jj); tau(jj,4+(ii-1)*numTauVariables) = nan; end;

                tau(jj,5+(ii-1)*numTauVariables) = mean(vPF_P.*vPF_P, 'omitmissing'); tauHeader{5+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :vPF'vPF'");
                if rotatedSonFlag(jj); tau(jj,5+(ii-1)*numTauVariables) = nan; end;

                tau(jj,6+(ii-1)*numTauVariables) = mean(wPF_P.*wPF_P, 'omitmissing'); tauHeader{6+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :wPF'wPF'");
                if rotatedSonFlag(jj); tau(jj,6+(ii-1)*numTauVariables) = nan; end;

                tau(jj,7+(ii-1)*numTauVariables) = mean(uPF_P.*vPF_P, 'omitmissing'); tauHeader{7+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :uPF'vPF'");
                if rotatedSonFlag(jj); tau(jj,7+(ii-1)*numTauVariables) = nan; end;

                tau(jj,8+(ii-1)*numTauVariables) = mean(uPF_P.*wPF_P, 'omitmissing'); tauHeader{8+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :uPF'wPF'");
                if rotatedSonFlag(jj); tau(jj,8+(ii-1)*numTauVariables) = nan; end;

                tau(jj,9+(ii-1)*numTauVariables) = mean(vPF_P.*wPF_P, 'omitmissing'); tauHeader{9+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :vPF'wPF'");
                if rotatedSonFlag(jj); tau(jj,9+(ii-1)*numTauVariables) = nan; end;

                tau(jj,10+(ii-1)*numTauVariables) = mean(u_tilt_P.*u_tilt_P, 'omitmissing'); tauHeader{10+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :uTilt'uTilt'");
                if rotatedSonFlag(jj); tau(jj,10+(ii-1)*numTauVariables) = nan; end;

                tau(jj,11+(ii-1)*numTauVariables) = mean(v_tilt_P.*v_tilt_P, 'omitmissing'); tauHeader{11+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :vTilt'vTilt'");
                if rotatedSonFlag(jj); tau(jj,11+(ii-1)*numTauVariables) = nan; end;

                tau(jj,12+(ii-1)*numTauVariables) = mean(w_tilt_P.*w_tilt_P, 'omitmissing'); tauHeader{12+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :wTilt'wTilt'");
                if rotatedSonFlag(jj); tau(jj,12+(ii-1)*numTauVariables) = nan; end;

                tau(jj,13+(ii-1)*numTauVariables) = mean(u_tilt_P.*v_tilt_P, 'omitmissing'); tauHeader{13+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :uTilt'vTilt'");
                if rotatedSonFlag(jj); tau(jj,13+(ii-1)*numTauVariables) = nan; end;

                tau(jj,14+(ii-1)*numTauVariables) = mean(u_tilt_P.*w_tilt_P, 'omitmissing'); tauHeader{14+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :uTilt'wTilt'");
                if rotatedSonFlag(jj); tau(jj,14+(ii-1)*numTauVariables) = nan; end;

                tau(jj,15+(ii-1)*numTauVariables) = mean(v_tilt_P.*w_tilt_P, 'omitmissing'); tauHeader{15+(ii-1)*numTauVariables} = strcat(num2str(sonHeight),"m :vTilt'wTilt'");
                if rotatedSonFlag(jj); tau(jj,15+(ii-1)*numTauVariables) = nan; end;
                
                % find tke
                tke(jj,2+ii) = 1/2*(mean(uP.^2, 'omitmissing')+mean(vP.^2, 'omitmissing')+mean(wP.^2, 'omitmissing')); tkeHeader{2+ii} = strcat(num2str(sonHeight),'m :0.5(u''^2+v''^2+w''^2)');
                if rotatedSonFlag(jj); tke(jj,2+ii) = nan; end;

                % find turbulent transport of TKE
                turbtr(jj,2+ii) = mean(wPF_P .* 1/2.*(uPF_P.^2 + vPF_P.^2+ wPF_P.^2), 'omitmissing'); turbtrHeader{2+ii} = strcat(num2str(sonHeight),"m :mean(w'e')");
                if rotatedSonFlag(jj); turbtr(jj,2+ii) = nan; end;
                
                if info.calcDissipation
                    % find dissipation rate of TKE
                    epsilon(jj,1+ii) = calc_dissipation_rate(uPF_P, mean(uPF(bp(jj)+1:bp(jj+1)), 'omitmissing'), 1./info.tableScanFrequency); epsilonHeader{1+ii} = strcat(num2str(sonHeight),"m :epsilon");
                    if rotatedSonFlag(jj); epsilon(jj,1+ii) = nan; end;
                end

                % find heat fluxes in SNSP system
                numH_SNSPVariables = 4;
                H_SNSP(jj,2+(ii-1)*numH_SNSPVariables) = mean(u_tilt_P.*Theta_v_sonP, 'omitmissing'); H_SNSPHeader{2+(ii-1)*numH_SNSPVariables} = strcat(num2str(sonHeight),'m :uTHv');
                if rotatedSonFlag(jj); H_SNSP(jj,2+(ii-1)*numH_SNSPVariables) = nan; end;
                
                H_SNSP(jj,3+(ii-1)*numH_SNSPVariables) = mean(v_tilt_P.*Theta_v_sonP, 'omitmissing'); H_SNSPHeader{3+(ii-1)*numH_SNSPVariables} = strcat(num2str(sonHeight),'m :vTHv');
                if rotatedSonFlag(jj); H_SNSP(jj,3+(ii-1)*numH_SNSPVariables) = nan; end;
                
                H_SNSP(jj,4+(ii-1)*numH_SNSPVariables) = mean(w_tilt_P.*Theta_v_sonP, 'omitmissing'); H_SNSPHeader{4+(ii-1)*numH_SNSPVariables} = strcat(num2str(sonHeight),'m :wTHv');
                if rotatedSonFlag(jj); H_SNSP(jj,4+(ii-1)*numH_SNSPVariables) = nan; end;

                % compute alpha_u and alpha_v for total vertical heat fluxes
                [alpha_u, alpha_v] = calc_SNSP_angle(direction(jj), angle);
                H_SNSP(jj,5+(ii-1)*numH_SNSPVariables) = mean(w_tilt_P.*Theta_v_sonP, 'omitmissing').*cosd(angle) - mean(u_tilt_P.*Theta_v_sonP, 'omitmissing').*sind(alpha_u) - mean(v_tilt_P.*Theta_v_sonP, 'omitmissing').*sind(alpha_v); H_SNSPHeader{5+(ii-1)*numH_SNSPVariables} = strcat(num2str(sonHeight),'m :wTHv_vert');
                if rotatedSonFlag(jj); H_SNSP(jj,5+(ii-1)*numH_SNSPVariables) = nan; end;

                
                % find Ts'w' and Ts'wPF' from sonic
                H(jj,4+(ii-1)*12) = mean(wP.*TsonP, 'omitmissing'); Hheader{4+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''w''');
                if unrotatedSonFlag(jj)||TsonFlag(jj); H(jj,4+(ii-1)*12) = nan; end;
                
                H(jj,5+(ii-1)*12) = mean(wPF_P.*Theta_v_sonP, 'omitmissing'); Hheader{5+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Theta_v''wPF'''); % Modified by Diane 2025/08/04 to account for slope angle impacts on vertical heat flux
                if rotatedSonFlag(jj)||TsonFlag(jj); H(jj,5+(ii-1)*12) = nan; end;
                
                Hlat(jj,2+(ii-1)*12) = mean(uP.*TsonP, 'omitmissing');HlatHeader{2+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''u''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,2+(ii-1)*12) = nan; end;
                
                Hlat(jj,3+(ii-1)*12) = mean(uPF_P.*TsonP, 'omitmissing');HlatHeader{3+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''uPF''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,3+(ii-1)*12) = nan; end;
                
                Hlat(jj,4+(ii-1)*12) = mean(vP.*TsonP, 'omitmissing');HlatHeader{4+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''v''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,4+(ii-1)*12) = nan; end;
                
                Hlat(jj,5+(ii-1)*12) = mean(vPF_P.*TsonP, 'omitmissing');HlatHeader{5+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Ts''vPF''');
                if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,5+(ii-1)*12) = nan; end;
                
                % find obukhov length
                kappa = 0.4;
                g = 9.81;
                T0_L = mean(thetaSon(bp(jj)+1:bp(jj+1)), "omitmissing") + 273.15; % Theta_v in K; Modified by Diane 2025/08/05
                % T0_L = Tref_Kavg(jj);
                uStarCubed = tau(jj,3+(ii-1)*numTauVariables).^(3/2); % sqrt(uPF'*wPF'+vPF'*wPF')
                % H0_L = H(jj,5+(ii-1)*12); % Theta_v''wPF''cos_alpha-Theta_v''uPF''sin_alpha
                H0_L = H_SNSP(jj,5+(ii-1)*numH_SNSPVariables); % use total vertical heat fluxes
                
                L(jj,2+(ii-1)) = -uStarCubed/(kappa*g/T0_L*H0_L); Lheader{2+(ii-1)} = strcat(num2str(sonHeight),'m L:sqrt(uPF''*wPF''+vPF''*wPF'')^3/2*Th_v/(k*g*wThv_vert)');
                if rotatedSonFlag(jj)||TsonFlag(jj); L(jj,2+(ii-1)) = nan; end;
                
                % store thetaSonP in raw structure
                if info.saveRawConditionedData
                    raw.sonTs(bp(jj)+1:bp(jj+1),ii) = Tson(bp(jj)+1:bp(jj+1));
%                     raw.sonTsPrime(bp(jj)+1:bp(jj+1),ii) = TsonP;
                    % raw.Theta_v_son(bp(jj)+1:bp(jj+1),ii) = thetaSonAir(bp(jj)+1:bp(jj+1));
                    raw.Theta_v_son(bp(jj)+1:bp(jj+1),ii) = thetaSon(bp(jj)+1:bp(jj+1)); % Modified by Diane to get sonic virtual potential temperature
%                     raw.Theta_v_sonPrime(bp(jj)+1:bp(jj+1),ii) = Theta_v_sonP;
                end
                
                % find Th'w' and Th'wPF' from sonic
                if ~isempty(thetaSonAir)
                    thetaSonAirP = nandetrend(thetaSonAir(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,6+(ii-1)*12) = mean(wP.*thetaSonAirP, 'omitmissing'); Hheader{6+(ii-1)*12} = strcat(num2str(sonHeight),'m son:T_air''w'''); % modified by Diane
                    if unrotatedSonFlag(jj)||TsonFlag(jj); H(jj,6+(ii-1)*12) = nan; end;
                    
                    H(jj,7+(ii-1)*12) = mean(wPF_P.*thetaSonAirP, 'omitmissing');Hheader{7+(ii-1)*12} = strcat(num2str(sonHeight),'m son:T_air''wPF'''); % modified by Diane
                    if rotatedSonFlag(jj)||TsonFlag(jj); H(jj,7+(ii-1)*12) = nan; end;
                end
                
                % find T'w' and T'wPF' from fw
                if ~isempty(fw)
                    fwP = nandetrend(fw(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,8+(ii-1)*12) = mean(wP.*fwP, 'omitmissing'); Hheader{8+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:T''w''');
                    if unrotatedSonFlag(jj)||fwFlag(jj); H(jj,8+(ii-1)*12) = nan; end;
                    
                    H(jj,9+(ii-1)*12) = mean(wPF_P.*fwP, 'omitmissing');Hheader{9+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:T''wPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); H(jj,9+(ii-1)*12) = nan; end;
                end
                
                % find Th_s'w' and Th_s'wPF' from sonics
                if ~isempty(thetaSon)
                    thetaSonP = nandetrend(thetaSon(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,10+(ii-1)*12) = mean(wP.*thetaSonP, 'omitmissing'); Hheader{10+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''w''');
                    if unrotatedSonFlag(jj)||TsonFlag(jj); H(jj,10+(ii-1)*12) = nan; end;
                    
                    H(jj,11+(ii-1)*12) = mean(wPF_P.*thetaSonP, 'omitmissing');Hheader{11+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''wPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); H(jj,11+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,6+(ii-1)*12) = mean(uP.*thetaSonP, 'omitmissing'); HlatHeader{6+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''u''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,6+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,7+(ii-1)*12) = mean(uPF_P.*thetaSonP, 'omitmissing');HlatHeader{7+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''uPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,7+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,8+(ii-1)*12) = mean(vP.*thetaSonP, 'omitmissing'); HlatHeader{8+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''v''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,8+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,9+(ii-1)*12) = mean(vPF_P.*thetaSonP, 'omitmissing');HlatHeader{9+(ii-1)*12} = strcat(num2str(sonHeight),'m son:Th_s''vPF''');
                    if rotatedSonFlag(jj)||TsonFlag(jj); Hlat(jj,9+(ii-1)*12) = nan; end;
                end
                
                % find Th'w' and Th'wPF' from fw
                if ~isempty(thetaFw)
                    thetaFwP = nandetrend(thetaFw(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,12+(ii-1)*12) = mean(wP.*thetaFwP, 'omitmissing'); Hheader{12+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''w''');
                    if unrotatedSonFlag(jj)||fwFlag(jj); H(jj,12+(ii-1)*12) = nan; end;
                    
                    H(jj,13+(ii-1)*12) = mean(wPF_P.*thetaFwP, 'omitmissing');Hheader{13+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''wPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); H(jj,13+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,10+(ii-1)*12) = mean(uP.*thetaFwP, 'omitmissing'); HlatHeader{10+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''u''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,10+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,11+(ii-1)*12) = mean(uPF_P.*thetaFwP, 'omitmissing');HlatHeader{11+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''uPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,11+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,12+(ii-1)*12) = mean(vP.*thetaFwP, 'omitmissing'); HlatHeader{12+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''v''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,12+(ii-1)*12) = nan; end;
                    
                    Hlat(jj,13+(ii-1)*12) = mean(vPF_P.*thetaFwP, 'omitmissing');HlatHeader{13+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:Th''vPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); Hlat(jj,13+(ii-1)*12) = nan; end;
                    
                    % store fwthetaP in raw structure
                    if info.saveRawConditionedData
                        raw.fwThPrime(bp(jj)+1:bp(jj+1),ii) = thetaFwP;
                        raw.fwTh(bp(jj)+1:bp(jj+1),ii) = thetaFw(bp(jj)+1:bp(jj+1));
                        raw.fwT(bp(jj)+1:bp(jj+1),ii) = fw(bp(jj)+1:bp(jj+1));
                    end
                end
                
                % find VTh'w' and VTh'wPF' from fw
                if ~isempty(VthetaFw)
                    VthetaFwP = nandetrend(VthetaFw(bp(jj)+1:bp(jj+1)),info.detrendingFormat);
                    
                    H(jj,14+(ii-1)*12) = mean(wP.*VthetaFwP, 'omitmissing'); Hheader{14+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:VTh''w''');
                    if unrotatedSonFlag(jj)||fwFlag(jj); H(jj,14+(ii-1)*12) = nan; end;
                    
                    H(jj,15+(ii-1)*12) = mean(wPF_P.*VthetaFwP, 'omitmissing');Hheader{15+(ii-1)*12} = strcat(num2str(sonHeight),'m fw:VTh''wPF''');
                    if rotatedSonFlag(jj)||fwFlag(jj); H(jj,15+(ii-1)*12) = nan; end;
                end
                
                % find latent heat flux if it exists at height
                if ~isempty(irgaH2O) || ~isempty(KH2O)
                    
                    % find H2O measurement
                    if ~isempty(irgaH2O)
                        rho_v = irgaH2O; % g/m^3
                        H2Otype = 0;
                        rho_CO2 = irgaCO2; % mg/m^3
                        if isfield(sensorInfo,'irgaH2O') 
                            H2OsensorNumber = find(abs(sensorInfo.irgaH2O(:,3)-sonHeight)<0.2);  % find IRGA number to store raw pertubations
                        end
                        if isfield(sensorInfo,'LiH2O') %&& isempty(H2OsensorNumber)
                            H2OsensorNumber = find(abs(sensorInfo.LiH2O(:,3)-sonHeight)<0.2);  % find IRGA number to store raw pertubations
                        end
                        CO2sensorNumber = H2OsensorNumber;
                    else
                        rho_v = KH2O; % g/m^3
                        H2Otype = 1;
                        H2OsensorNumber = find(abs(sensorInfo.KH2O(:,3)-sonHeight)<0.2);  % find KH2O number to store raw pertubations
                    end
                    
                    % declare constants
                    Mv = 18.0153;   % Molar Mass of H2O (g/mol)
                    Md = 28.97;     % Molar Mass of Dry Air (g/mol)
                    
                    % find heat flux for WPL (wPF_P*TsonP)
                    kinSenFlux = H(jj,5+(ii-1)*12); % K m/s
                    
                    % find the latent heat of vaporizatoin
                    Lv = (2.501-0.00237*(Tref_Kavg(jj)-273.15))*10^3; % Latent Heat of Vaporization (J/g) Stoll P. 641
                    
                    % find H2O pertubations in g/m^3
                    H2Op = nandetrend(rho_v(bp(jj)+1:bp(jj+1)),info.detrendingFormat); %(g/m^3)

                    % More sigma computation
                    sigma(jj,11+(ii-1)*numSigmaVariables) = std(rho_v(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{11+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_H2O');
                    if TsonFlag(jj) || rotatedSonFlag(jj); sigma(jj,11+(ii-1)*numSigmaVariables) = nan; end;
                
                    
                    % more flux correlation
                    fuw = uPF_P.*wPF_P; fHw = H2Op.*wPF_P; indx = ~(isnan(fuw) | isnan(fHw));
                    R_mat_temp = corrcoef(fuw(indx), fHw(indx));
                    R(jj, 3+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{3+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_uPFwPF_wPFH2O');
                    if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 3+(ii-1)*numRVariables) = nan; end;
                    
                    fHw = H2Op.*wPF_P; fwt = Theta_v_sonP.*wPF_P; indx = ~(isnan(fHw) | isnan(fwt));
                    R_mat_temp = corrcoef(fHw(indx), fwt(indx));
                    R(jj, 4+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{4+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPFH2O_wPFThetav');
                    if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 4+(ii-1)*numRVariables) = nan; end;

                    % more transport efficiency
                    eta(jj, 4+(ii-1)*numetaVariables) = find_eta(wPF_P, H2Op); etaHeader{4+(ii-1)*numetaVariables} = strcat(num2str(sonHeight),'m :eta_wPFH2O');
                    if TsonFlag(jj) || rotatedSonFlag(jj); eta(jj, 4+(ii-1)*numetaVariables) = nan; end;
                    
                    delta_flux_ctrb(jj, 4+(ii-1)*numDelta_flux_ctrbVariables) = find_delta_flux(wPF_P, H2Op); delta_flux_ctrbHeader{4+(ii-1)*numDelta_flux_ctrbVariables} = strcat(num2str(sonHeight),'m :S_wPFH2O');
                    if TsonFlag(jj) || rotatedSonFlag(jj); delta_flux_ctrb(jj, 4+(ii-1)*numDelta_flux_ctrbVariables) = nan; end;
                    
                    delta_time_ctrb(jj, 4+(ii-1)*numDelta_time_ctrbVariables) = find_delta_time(wPF_P, H2Op); delta_time_ctrbHeader{4+(ii-1)*numDelta_time_ctrbVariables} = strcat(num2str(sonHeight),'m :D_wPFH2O');
                    if TsonFlag(jj) || rotatedSonFlag(jj); delta_time_ctrb(jj, 4+(ii-1)*numDelta_time_ctrbVariables) = nan; end;
                

                    % store H2OP in raw structure
                    if info.saveRawConditionedData
                        raw.rhovPrime(bp(jj)+1:bp(jj+1),H2OsensorNumber) = H2Op; % (g/m^3)
                        raw.rhov(bp(jj)+1:bp(jj+1),H2OsensorNumber) = rho_v(bp(jj)+1:bp(jj+1)); % (g/m^3)
                    end
                    
                    E = mean(wP.*H2Op, 'omitmissing');     % [g /m^2/s] find evaporation flux from unrotated data
                    EPF = mean(wPF_P.*H2Op, 'omitmissing');  % [g /m^2/s] find evaporation flux from rotated data
                    
                    LHflux(jj,2+(ii-1)*7) = Lv; LHfluxHeader{2+(ii-1)*7} = strcat(num2str(sonHeight),'m Lv(J/g)');
                    
                    LHflux(jj,3+(ii-1)*7) = E; LHfluxHeader{3+(ii-1)*7} = strcat(num2str(sonHeight),'m w'':E(g/m^2s)');
                    if unrotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,3+(ii-1)*7) = nan; end;
                    
                    LHflux(jj,4+(ii-1)*7) = EPF; LHfluxHeader{4+(ii-1)*7} = strcat(num2str(sonHeight),'m wPF'':E(g/m^2s)');
                    if rotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,4+(ii-1)*7) = nan; end;
                    
                    % apply O2 correction to KH2O
                    if H2Otype % see: https://www.eol.ucar.edu/content/corrections-sensible-and-latent-heat-flux-measurements 
                               % EVERYTHING ELSE IS DEFINITELY WRONG!!! THIS LOOKS RIGHT
                        ko = -0.0045; % effective O2 absorption coefficient.  Value from: Tanner et al., 1993
                        kw = -0.153;  % effective H2O absorption coefficient.  Value from instrument calibration.  -0.153 is from 2012 calibration of the Pardyjak lab KH2O
                        CkO = 0.23*ko/kw;
                        O2correction = CkO*rhodAvg(jj)/Tref_Kavg(jj)*kinSenFlux*1000;  % multiply by 1000 to put in g/m^2/s
                        E = E+O2correction;
                        EPF = EPF+O2correction;
                        LHflux(jj,7+(ii-1)*7) = Lv*E; LHfluxHeader{7+(ii-1)*7} = strcat(num2str(sonHeight),'m O2 no WPL,w'' (W/m^2)');
                        LHflux(jj,8+(ii-1)*7) = Lv*EPF; LHfluxHeader{8+(ii-1)*7} = strcat(num2str(sonHeight),'m O2 no WPL,wPF'' (W/m^2)');
                        
                        % create WPL headers to include O2 correction
                        LHfluxHeader{5+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL and O2, w'' (W/m^2)');
                        LHfluxHeader{6+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL and O2, wPF'' (W/m^2)');
                    end
                    
                    % apply WPL Corrections (E.C. by Marc Aubinet 97).  Headers are created above for KH2O and below for EC150
                    LHflux(jj,5+(ii-1)*7) = 1000*Lv*(1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(E./1000+(rhovAvg(jj)/Tref_Kavg(jj))*kinSenFlux); 
                    if unrotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,5+(ii-1)*7) = nan; end;
                    
                    LHflux(jj,6+(ii-1)*7) = 1000*Lv*(1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(EPF./1000+(rhovAvg(jj)/Tref_Kavg(jj))*kinSenFlux); 
                    if unrotatedSonFlag(jj)||H2OFlag(jj); LHflux(jj,6+(ii-1)*7) = nan; end;
                    
                    
                    if ~H2Otype  % find CO2 flux and create H2O WPL headers
                        
                        % create WPL, H2O headers 
                        LHfluxHeader{5+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL, w'' (W/m^2)');
                        LHfluxHeader{6+(ii-1)*7} = strcat(num2str(sonHeight),'m WPL, wPF'' (W/m^2)');
                        
                        rho_CO2p = nandetrend(rho_CO2(bp(jj)+1:bp(jj+1)),info.detrendingFormat)./1e6; % kg/m^3
                        rho_CO2avg = mean(rho_CO2(bp(jj)+1:bp(jj+1)), 'omitmissing')/1e6; % kg/m^3
                        evapFlux = LHflux(jj,6+(ii-1)*7)/Lv/1000; % WPL, wPF'' (kg/m^2s)

                        % compute external CO2 fluctuations induced by water vapor and temperature
                        rho_externalp = Md./Mv.*(rho_CO2avg./rhodAvg(jj)).*(H2Op./1000) + rho_CO2avg.*(1+Md./Mv .* rhovAvg(jj)./rhodAvg(jj)).*T_air_sonP./(T_air_sonM+273.15); % following Detto and Katul 2007, Code by Diane

                        CO2flux(jj,2+(ii-1)*4) = mean(wP.*rho_CO2p, 'omitmissing'); CO2fluxHeader{2+(ii-1)*4} = strcat(num2str(sonHeight),'m w'':CO2(kg/m^2s)');
                        if unrotatedSonFlag(jj)||CO2Flag(jj); CO2flux(jj,2+(ii-1)*4) = nan; end;
                        
                        CO2flux(jj,3+(ii-1)*4) = mean(wPF_P.*rho_CO2p, 'omitmissing'); CO2fluxHeader{3+(ii-1)*4} = strcat(num2str(sonHeight),'m wPF'':CO2(kg/m^2s)');
                        if rotatedSonFlag(jj)||CO2Flag(jj); CO2flux(jj,3+(ii-1)*4) = nan; end;
                        
                        CO2flux(jj,4+(ii-1)*4) = CO2flux(jj,2+(ii-1)*4) + Md/Mv*(rho_CO2avg/rhodAvg(jj))*evapFlux + (1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(rho_CO2avg/Tref_Kavg(jj))*kinSenFlux; CO2fluxHeader{4+(ii-1)*4} = strcat(num2str(sonHeight),'m WPL,w'':CO2(kg/m^2s)');
                        
                        % CO2flux(jj,5+(ii-1)*4) = CO2flux(jj,3+(ii-1)*4) + Md/Mv*(rho_CO2avg/rhodAvg(jj))*evapFlux + (1+Md/Mv*rhovAvg(jj)/rhodAvg(jj))*(rho_CO2avg/Tref_Kavg(jj))*kinSenFlux; CO2fluxHeader{5+(ii-1)*4} = strcat(num2str(sonHeight),'m WPL,wPF'':CO2(mol/m^2s)');
                        CO2flux(jj,5+(ii-1)*4) = CO2flux(jj,3+(ii-1)*4) + mean(wPF_P.*rho_externalp, 'omitmissing'); CO2fluxHeader{5+(ii-1)*4} = strcat(num2str(sonHeight),'m: wPF''CO2_WPL''(m/s kg/m^3)'); % Modified by Diane
                        
                        % More sigma computation
                        sigma(jj,12+(ii-1)*numSigmaVariables) = std(rho_CO2(bp(jj)+1:bp(jj+1)), 'omitmissing'); sigmaHeader{12+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_CO2');
                        if TsonFlag(jj) || rotatedSonFlag(jj); sigma(jj,12+(ii-1)*numSigmaVariables) = nan; end;

                        sigma(jj,13+(ii-1)*numSigmaVariables) = std(rho_CO2(bp(jj)+1:bp(jj+1))+rho_externalp, 'omitmissing'); sigmaHeader{13+(ii-1)*numSigmaVariables} = strcat(num2str(sonHeight),'m :sigma_CO2_WPL');
                        if TsonFlag(jj) || rotatedSonFlag(jj); sigma(jj,13+(ii-1)*numSigmaVariables) = nan; end;

                
                        % store CO2P in raw structure
                        if info.saveRawConditionedData
                            raw.rhoCO2Prime(bp(jj)+1:bp(jj+1),CO2sensorNumber) = rho_CO2p;
                            raw.rhoCO2exteralPrime(bp(jj)+1:bp(jj+1),CO2sensorNumber) = rho_externalp;
                            raw.rhoCO2(bp(jj)+1:bp(jj+1),CO2sensorNumber) = rho_CO2(bp(jj)+1:bp(jj+1));
                        end

                        % more flux correlation
                        fuw = uPF_P.*wPF_P; fCw = (rho_CO2p + rho_externalp).*wPF_P; indx = ~(isnan(fuw) | isnan(fCw));
                        R_mat_temp = corrcoef(fuw(indx), fCw(indx));
                        R(jj, 5+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{5+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_uPFwPF_wPFCO2_WPL');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 5+(ii-1)*numRVariables) = nan; end;
                    
                        fCw = (rho_CO2p + rho_externalp).*wPF_P; fwt = Theta_v_sonP.*wPF_P; indx = ~(isnan(fCw) | isnan(fwt));
                        R_mat_temp = corrcoef(fCw(indx), fwt(indx));
                        R(jj, 6+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{6+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPFCO2_WPL_wPFThetav');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 6+(ii-1)*numRVariables) = nan; end;
                        
                        indx = ~(isnan(wPF_P) | isnan(uPF_P));
                        R_mat_temp = corrcoef(wPF_P(indx), uPF_P(indx));
                        R(jj, 7+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{7+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPF_uPF');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 7+(ii-1)*numRVariables) = nan; end;
                        
                        indx = ~(isnan(wPF_P) | isnan(Theta_v_sonP));
                        R_mat_temp = corrcoef(wPF_P(indx), Theta_v_sonP(indx));
                        R(jj, 8+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{8+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPF_Theta_v');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 8+(ii-1)*numRVariables) = nan; end;
                        
                        indx = ~(isnan(wPF_P) | isnan(H2Op));
                        R_mat_temp = corrcoef(wPF_P(indx), H2Op(indx));
                        R(jj, 9+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{9+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPF_H2O');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 9+(ii-1)*numRVariables) = nan; end;
                        
                        indx = ~(isnan(wPF_P) | isnan(rho_CO2p));
                        R_mat_temp = corrcoef(wPF_P(indx), rho_CO2p(indx));
                        R(jj, 10+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{10+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPF_CO2');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 10+(ii-1)*numRVariables) = nan; end;
                        
                        indx = ~(isnan(wPF_P) | isnan(rho_CO2p+rho_externalp));
                        R_mat_temp = corrcoef(wPF_P(indx), rho_CO2p(indx)+rho_externalp(indx));
                        R(jj, 11+(ii-1)*numRVariables) = R_mat_temp(1,2); RHeader{11+(ii-1)*numRVariables} = strcat(num2str(sonHeight),'m :R_wPF_CO2_WPL');
                        if TsonFlag(jj) || rotatedSonFlag(jj); R(jj, 11+(ii-1)*numRVariables) = nan; end;
                        
                        % more transport efficiency
                        eta(jj, 5+(ii-1)*numetaVariables) = find_eta(wPF_P, rho_CO2p+rho_externalp); etaHeader{5+(ii-1)*numetaVariables} = strcat(num2str(sonHeight),'m :eta_wPFCO2_WPL');
                        if TsonFlag(jj) || rotatedSonFlag(jj); eta(jj, 5+(ii-1)*numetaVariables) = nan; end;
                        
                        delta_flux_ctrb(jj, 5+(ii-1)*numDelta_flux_ctrbVariables) = find_delta_flux(wPF_P, rho_CO2p+rho_externalp); delta_flux_ctrbHeader{5+(ii-1)*numDelta_flux_ctrbVariables} = strcat(num2str(sonHeight),'m :S_wPFCO2_WPL');
                        if TsonFlag(jj) || rotatedSonFlag(jj); delta_flux_ctrb(jj, 5+(ii-1)*numDelta_flux_ctrbVariables) = nan; end;
                        
                        delta_time_ctrb(jj, 5+(ii-1)*numDelta_time_ctrbVariables) = find_delta_time(wPF_P, rho_CO2p+rho_externalp); delta_time_ctrbHeader{5+(ii-1)*numDelta_time_ctrbVariables} = strcat(num2str(sonHeight),'m :D_wPFCO2_WPL');
                        if TsonFlag(jj) || rotatedSonFlag(jj); delta_time_ctrb(jj, 5+(ii-1)*numDelta_time_ctrbVariables) = nan; end;
                        
                        % compute skewness
                        numSkewVariables = 7;
                        skew(jj,2+(ii-1)*numSkewVariables) = skewness(uPF_P,0); skewHeader{2+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_uPF');
                        if rotatedSonFlag(jj); skew(jj,2+(ii-1)*numSkewVariables) = nan; end;
                        
                        skew(jj,3+(ii-1)*numSkewVariables) = skewness(vPF_P,0); skewHeader{3+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_vPF');
                        if rotatedSonFlag(jj); skew(jj,3+(ii-1)*numSkewVariables) = nan; end;
                        
                        skew(jj,4+(ii-1)*numSkewVariables) = skewness(wPF_P,0); skewHeader{4+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_wPF');
                        if rotatedSonFlag(jj); skew(jj,4+(ii-1)*numSkewVariables) = nan; end;
                        
                        skew(jj,5+(ii-1)*numSkewVariables) = skewness(Theta_v_sonP,0); skewHeader{5+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_Theata_v');
                        if rotatedSonFlag(jj); skew(jj,5+(ii-1)*numSkewVariables) = nan; end;
                        
                        skew(jj,6+(ii-1)*numSkewVariables) = skewness(H2Op,0); skewHeader{6+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_H2O');
                        if rotatedSonFlag(jj); skew(jj,6+(ii-1)*numSkewVariables) = nan; end;
                        
                        skew(jj,7+(ii-1)*numSkewVariables) = skewness(rho_CO2p,0); skewHeader{7+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_CO2');
                        if rotatedSonFlag(jj); skew(jj,7+(ii-1)*numSkewVariables) = nan; end;
                        
                        skew(jj,8+(ii-1)*numSkewVariables) = skewness(rho_CO2p+rho_externalp,0); skewHeader{8+(ii-1)*numSkewVariables} = strcat(num2str(sonHeight),'m :skew_CO2_WPL');
                        if rotatedSonFlag(jj); skew(jj,8+(ii-1)*numSkewVariables) = nan; end;
                        
                        % compute lateral scalar fluxes
                        numFLatVariables = 4;
                        Flux_lat(jj,2+(ii-1)*numFLatVariables) = mean(uPF_P.*Theta_v_sonP, 'omitmissing'); Flux_latHeader{2+(ii-1)*numFLatVariables} = strcat(num2str(sonHeight),"m :uPF'theta_v'");
                        if rotatedSonFlag(jj); Flux_lat(jj,2+(ii-1)*numFLatVariables) = nan; end;
                        
                        Flux_lat(jj,3+(ii-1)*numFLatVariables) = mean(uPF_P.*H2Op, 'omitmissing'); Flux_latHeader{3+(ii-1)*numFLatVariables} = strcat(num2str(sonHeight),"m :uPF'H2O'");
                        if rotatedSonFlag(jj); Flux_lat(jj,3+(ii-1)*numFLatVariables) = nan; end;
                        
                        Flux_lat(jj,4+(ii-1)*numFLatVariables) = mean(uPF_P.*rho_CO2p, 'omitmissing'); Flux_latHeader{4+(ii-1)*numFLatVariables} = strcat(num2str(sonHeight),"m :uPF'CO2'");
                        if rotatedSonFlag(jj); Flux_lat(jj,4+(ii-1)*numFLatVariables) = nan; end;
                        
                        Flux_lat(jj,5+(ii-1)*numFLatVariables) = mean(uPF_P.*(rho_CO2p+rho_externalp), 'omitmissing'); Flux_latHeader{5+(ii-1)*numFLatVariables} = strcat(num2str(sonHeight),"m :uPF'CO2_WPL'");
                        if rotatedSonFlag(jj); Flux_lat(jj,5+(ii-1)*numFLatVariables) = nan; end;
                        
                    end
                end
                
                
            end
        catch err
            message = strcat(err.message,'@ line',num2str(err.stack.line),' Problem with sonic at ',num2str(sonHeight),' m will be skipped');
            warning(message)
            if isempty(output.warnings{1})
                output.warnings{1,1} = message;
            else
                output.warnings{end+1,1} = message;
            end
        end
    end
    %------------- STORE OUTPUTS
    flag = logical(any(H,1)+isnan(H(1,:)));
    output.H = H(:,flag);
    output.Hheader = Hheader(flag);
    
    flag = logical(any(Hlat,1)+isnan(Hlat(1,:)));
    output.Hlat = Hlat(:,flag);
    output.HlatHeader = HlatHeader(flag);
    
    flag = logical(any(tau,1)+isnan(tau(1,:)));
    output.tau = tau(:,flag);
    output.tauHeader = tauHeader(flag);
    
    flag = logical(any(tke,1)+isnan(tke(1,:)));
    output.tke = tke(:,flag);
    output.tkeHeader = tkeHeader(flag);
    
    flag = logical(any(sigma,1)+isnan(sigma(1,:)));
    output.sigma = sigma(:,flag);
    output.sigmaHeader = sigmaHeader(flag);
    
    if info.storeExtraStats
        flag = logical(any(R,1)+isnan(R(1,:)));
        output.R = R(:,flag);
        output.RHeader = RHeader(flag);
        
        flag = logical(any(L,1)+isnan(L(1,:)));
        output.L = L(:,flag);
        output.Lheader = Lheader(flag);
    
        flag = logical(any(eta,1)+isnan(eta(1,:)));
        output.eta = eta(:,flag);
        output.etaHeader = etaHeader(flag);
    
        flag = logical(any(delta_flux_ctrb,1)+isnan(delta_flux_ctrb(1,:)));
        output.delta_flux_ctrb = delta_flux_ctrb(:,flag);
        output.delta_flux_ctrbHeader = delta_flux_ctrbHeader(flag);
    
        flag = logical(any(delta_time_ctrb,1)+isnan(delta_time_ctrb(1,:)));
        output.delta_time_ctrb = delta_time_ctrb(:,flag);
        output.delta_time_ctrbHeader = delta_time_ctrbHeader(flag);
        
        flag = logical(any(turbtr,1)+isnan(turbtr(1,:)));
        output.turbtr = turbtr(:,flag);
        output.turbtrHeader = turbtrHeader(flag);
    
        flag = logical(any(epsilon,1)+isnan(epsilon(1,:)));
        output.epsilon = epsilon(:,flag);
        output.epsilonHeader = epsilonHeader(flag);

        flag = logical(any(skew,1)+isnan(skew(1,:)));
        output.skew = skew(:,flag);
        output.skewHeader = skewHeader(flag);

        flag = logical(any(H_SNSP,1)+isnan(H_SNSP(1,:)));
        output.H_SNSP = H_SNSP(:,flag);
        output.H_SNSPHeader = H_SNSPHeader(flag);
    end

    flag = logical(any(Flux_lat,1)+isnan(Flux_lat(1,:)));
    output.Flux_lat = Flux_lat(:,flag);
    output.Flux_latHeader = Flux_latHeader(flag);
    
    if size(derivedT,2) > 1
        flag = logical(any(derivedT,1)+isnan(derivedT(1,:)));
        output.derivedT = derivedT(:,flag);
        output.derivedTheader = derivedTheader(flag);
    end
    
    if size(LHflux,2) > 1
        flag = logical(any(LHflux,1)+isnan(LHflux(1,:)));
        output.LHflux = LHflux(:,flag);
        output.LHfluxHeader = LHfluxHeader(flag);
    end
    
    if size(CO2flux,2) > 1
        flag = logical(any(CO2flux,1)+isnan(CO2flux(1,:)));
        output.CO2flux = CO2flux(:,flag);
        output.CO2fluxHeader = CO2fluxHeader(flag);
    end
catch err
    message = strcat(err.message,'@ line ',num2str(err.stack.line),' UNABLE TO FIND FLUXES AT All HEIGHTS');
    warning(message)
    raw = [];
    if isempty(output.warnings{1})
        output.warnings{1,1} = message;
    else
        output.warnings{end+1,1} = message;
    end
end
end