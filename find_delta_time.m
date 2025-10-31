function delta_time = find_delta_time(ww, uu)
    flux_all_events = ww .* uu;
    flux_total = sum(flux_all_events, 'omitmissing');
    flux_downgradient_flag = flux_all_events .* flux_total > 0;
    flux_ejection_flag = flux_downgradient_flag & ww>0;
    flux_sweep_flag = flux_downgradient_flag & ww<0;
    delta_time = sum(flux_ejection_flag) ./ length(ww) - sum(flux_sweep_flag) ./ length(ww);
end