function [virtualTheta, r, rho_airmoist, rhod, rhov] = get_virtulPotTemp(altitude, level, Tair, RH, P_air, usePelevation)
    % altitude: station altitude, @ 1980 m for Dolly Tower at FM site
    % level: HMP height
    
    if median(Tair, 'omitmissing')<200
        Tair = Tair + 273.15;  % put Tair in K
    end
    
    if usePelevation
        P_air = 101325*(1-2.25577*10^-5*(altitude+level))^5.25588; % Pa, find pRef from elevation: http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
        
    else
        if median(P_air, 'omitmissing')<1000
            P_air = P_air.*1000; % put P_air in Pa
        end
    end
    e_saturated = 1000 .* exp(52.57633 - 6790.4985./Tair - 5.02808.*log(Tair)); % Pa, find saturated water vapor presure from ATM 133 lecture
    e_air = e_saturated .* (RH./100); % Pa, unsaturated water vapor presure
    r = 621.97 .* e_air ./ (P_air - e_air); % g/kg, unsaturated mixing ratio

    % find potential temperature
    Gamma = 0.0098; % Dry Lapse Rate. K/m
    theta = Tair + Gamma .* level; % use height above ground level here
    % P0 = 101325;
    % theta = Tair .* ((P0 ./ P_air).^0.286);

    % Calculate virtual potential temperature
    virtualTheta = theta .* (1 + (0.61 .* r./1000)); % K

    % calculate moist air density
    % constants
    Rd = 287.058;  % [J/K/kg] Gas constant for air
    Rv = 461.495;  % [J/K/kg] Gas constant for water vapor
    rhod = (P_air - e_air)./(Rd*Tair); % density of dry air (kg/m^3)
    rhov = e_air./(Rv*Tair); % density of water vapor (kg/m^3)
    rho_airmoist = rhod+rhov; % density of moist air (kg/m^3)

end
