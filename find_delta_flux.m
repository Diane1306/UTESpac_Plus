function delta_flux = find_delta_flux(ww, uu)
    flux_all_events = ww .* uu;
    flux_total = sum(flux_all_events, 'omitmissing');
    flux_downgradient_flag = flux_all_events .* flux_total > 0;
    flux_ejection_flag = flux_downgradient_flag & ww>0;
    flux_sweep_flag = flux_downgradient_flag & ww<0;
    delta_flux = sum(flux_all_events(flux_ejection_flag), 'omitmissing') ./ flux_total - sum(flux_all_events(flux_sweep_flag), 'omitmissing') ./ flux_total;
end