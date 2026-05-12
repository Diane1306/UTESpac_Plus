function [a1, a2] = calc_SNSP_angle(phi,alpha)
    % phi is the meteorological wind direction from north in dgrees
    % alpha is the slope angle
    a3 = alpha;
    kesi = phi - 30; % French Meadows downslope wind is 30° from north
    a1 = asind(cosd(kesi).*sind(a3));
    a2 = asind(cosd(kesi-90).*sind(a3));
end