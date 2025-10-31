function eta = find_eta(ww, uu)
    flux_all_events = ww .* uu;
    flux_total = sum(flux_all_events, 'omitmissing');
    flux_downgradient_flag = flux_all_events .* flux_total > 0;
    eta = flux_total ./ sum(flux_all_events(flux_downgradient_flag), 'omitmissing');
end