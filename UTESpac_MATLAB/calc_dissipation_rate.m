function epsilon = calc_dissipation_rate(uP, Um, dt)
    % U is the mean velocity to convert time to distance
    y1 = calc_structure_function(uP);
    nPoints = length(uP);
    timemax_insec = nPoints*dt/2;
    r_value = linspace(dt.*Um, timemax_insec.*Um, floor(nPoints/2));
    y2 = r_value.^(2/3);
    [~, I] = min(abs(log(y2(1:20))-log(y1(1:20))));
    c_A2 = (y2 ./ (y2(I) ./ (y1(I)))) ./ y2;
    epsilon_array = (c_A2 .^ (3/2)) .* 0.35;
    epsilon = epsilon_array(1);
    
    function structure_function = calc_structure_function(u)
        % u is the velocity time series
        n = length(u);
        r = 1:floor(n/2);
        structure_function = zeros(size(r));
        for i = 1:length(r)
            ri = r(i);
            vel_diff = (u(1+ri:end) - u(1:end-ri)).^2;
            structure_function(i) = mean(vel_diff(:), 'omitmissing'); % average over the grid
        end
    end
end

