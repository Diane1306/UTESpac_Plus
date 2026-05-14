function [fluxQC_row, fluxQCHeader] = calc_SSITC_flags_ForestComplexTerrain( ...
    fluxQC_row, fluxQCHeader, ...
    jj, ii, sonHeight, info, ...
    wPF_P, uPF_P, vPF_P, Theta_v_sonP, ...
    H2Op, rhov_externalp, rho_CO2p, rhoc_externalp, ...
    tau, L, ...
    rotatedSonFlag, TsonFlag, H2OFlag, CO2Flag)

% calc_SSITC_flags_ForestComplexTerrain
%
% Outputs 0-1-2 quality flags for:
%   TAU_SSITC_TEST
%   TAU_SS_ONLY_TEST
%   H_SSITC_TEST
%   H_SS_ONLY_TEST
%   LE_SSITC_TEST
%   LE_SS_ONLY_TEST
%   FC_SSITC_TEST
%   FC_SS_ONLY_TEST
%
% SSITC flags:
%   steady-state test + canopy-aware w-based ITC test
%
% SS_ONLY flags:
%   steady-state test only, intended for slope-flow / complex-terrain
%   analyses where ITC assumptions are knowingly violated.
%
% Flag meaning:
%   0 = high quality
%   1 = moderate / usable
%   2 = poor / exclude or gap-fill

    numFluxQCVariables = 8;

    % ------------------------------------------------------------
    % Settings
    % ------------------------------------------------------------
    SSITC_subAvgMin = get_info_field(info, 'SSITC_subAvgMin', 5);
    canopyHeight    = get_info_field(info, 'canopyHeight', nan);
    displacementH   = get_info_field(info, 'displacementHeight', 0);
    useCanopyITC    = get_info_field(info, 'useCanopyITC', true);

    if mod(info.avgPer, SSITC_subAvgMin) ~= 0
        error('info.avgPer must be divisible by info.SSITC_subAvgMin.');
    end

    nSub = info.avgPer / SSITC_subAvgMin;

    % ------------------------------------------------------------
    % Compute ustar and z/L for w-based ITC
    % ------------------------------------------------------------
    numTauVariables = 14;
    tau_col = 3 + (ii-1)*numTauVariables;

    if size(tau,2) >= tau_col && ~isnan(tau(jj,tau_col)) && tau(jj,tau_col) >= 0
        ustar = sqrt(tau(jj,tau_col));
    else
        ustar = nan;
    end

    L_col = 2 + (ii-1);

    if size(L,2) >= L_col && ~isnan(L(jj,L_col)) && L(jj,L_col) ~= 0
        zeta = max(sonHeight - displacementH, 0.1) ./ L(jj,L_col);
    else
        zeta = nan;
    end

    sigma_w = std(wPF_P, 'omitmissing');

    itcDev_w = local_calc_ITCw_deviation( ...
        sigma_w, ustar, zeta, sonHeight, canopyHeight, useCanopyITC);

    % Column offset for this sonic
    base = 2 + (ii-1)*numFluxQCVariables;

    % ============================================================
    % 1. TAU flags
    % ============================================================
    tauSSDev_uw = local_calc_SS_deviation_full(uPF_P, wPF_P, nSub);
    tauSSDev_vw = local_calc_SS_deviation_full(vPF_P, wPF_P, nSub);
    tauSSDev = max([tauSSDev_uw, tauSSDev_vw], [], 'omitmissing');

    tauSSITCFlag = local_dev_to_SSITC_flag(tauSSDev, itcDev_w);
    tauSSonlyFlag = local_dev_to_SSonly_flag(tauSSDev);

    if rotatedSonFlag(jj)
        tauSSITCFlag = 2;
        tauSSonlyFlag = 2;
    end

    fluxQC_row(base)   = tauSSITCFlag;
    fluxQC_row(base+1) = tauSSonlyFlag;

    fluxQCHeader{base}   = sprintf('%gm:TAU_SSITC_TEST', sonHeight);
    fluxQCHeader{base+1} = sprintf('%gm:TAU_SS_ONLY_TEST', sonHeight);

    % ============================================================
    % 2. H flags
    % ============================================================
    HSSDev = local_calc_SS_deviation_full(wPF_P, Theta_v_sonP, nSub);

    HSSITCFlag = local_dev_to_SSITC_flag(HSSDev, itcDev_w);
    HSSonlyFlag = local_dev_to_SSonly_flag(HSSDev);

    if rotatedSonFlag(jj) || TsonFlag(jj)
        HSSITCFlag = 2;
        HSSonlyFlag = 2;
    end

    fluxQC_row(base+2) = HSSITCFlag;
    fluxQC_row(base+3) = HSSonlyFlag;

    fluxQCHeader{base+2} = sprintf('%gm:H_SSITC_TEST', sonHeight);
    fluxQCHeader{base+3} = sprintf('%gm:H_SS_ONLY_TEST', sonHeight);

    % ============================================================
    % 3. LE flags
    % ============================================================
    if ~isempty(H2Op) && ~isempty(rhov_externalp)
        H2O_WPL_P = H2Op + rhov_externalp .* 1e3; % g/m3 equivalent perturbation

        LESSDev = local_calc_SS_deviation_full(wPF_P, H2O_WPL_P, nSub);

        LESSITCFlag = local_dev_to_SSITC_flag(LESSDev, itcDev_w);
        LESSonlyFlag = local_dev_to_SSonly_flag(LESSDev);

        if rotatedSonFlag(jj) || H2OFlag(jj)
            LESSITCFlag = 2;
            LESSonlyFlag = 2;
        end
    else
        LESSITCFlag = nan;
        LESSonlyFlag = nan;
    end

    fluxQC_row(base+4) = LESSITCFlag;
    fluxQC_row(base+5) = LESSonlyFlag;

    fluxQCHeader{base+4} = sprintf('%gm:LE_SSITC_TEST', sonHeight);
    fluxQCHeader{base+5} = sprintf('%gm:LE_SS_ONLY_TEST', sonHeight);

    % ============================================================
    % 4. FC flags
    % ============================================================
    if ~isempty(rho_CO2p) && ~isempty(rhoc_externalp)
        CO2_WPL_P = rho_CO2p + rhoc_externalp; % kg/m3 equivalent perturbation

        FCSSDev = local_calc_SS_deviation_full(wPF_P, CO2_WPL_P, nSub);

        FCSSITCFlag = local_dev_to_SSITC_flag(FCSSDev, itcDev_w);
        FCSSonlyFlag = local_dev_to_SSonly_flag(FCSSDev);

        if rotatedSonFlag(jj) || CO2Flag(jj)
            FCSSITCFlag = 2;
            FCSSonlyFlag = 2;
        end
    else
        FCSSITCFlag = nan;
        FCSSonlyFlag = nan;
    end

    fluxQC_row(base+6) = FCSSITCFlag;
    fluxQC_row(base+7) = FCSSonlyFlag;

    fluxQCHeader{base+6} = sprintf('%gm:FC_SSITC_TEST', sonHeight);
    fluxQCHeader{base+7} = sprintf('%gm:FC_SS_ONLY_TEST', sonHeight);

end


% =====================================================================
% Helper functions
% =====================================================================

function value = get_info_field(info, fieldName, defaultValue)

    if isfield(info, fieldName)
        value = info.(fieldName);
    else
        value = defaultValue;
    end

end


function ssDev = local_calc_SS_deviation_full(x, w, nSub)

% Foken-style steady-state test.
% Compares full-period covariance with the mean of subperiod covariances.
% Subperiod means are removed within each subperiod.

    ssDev = nan;

    if isempty(x) || isempty(w) || numel(x) ~= numel(w)
        return
    end

    valid = ~(isnan(x) | isnan(w));

    if sum(valid) < 10
        return
    end

    % Full-period covariance
    xFull = x - mean(x, 'omitmissing');
    wFull = w - mean(w, 'omitmissing');

    covFull = mean(xFull .* wFull, 'omitmissing');

    if isnan(covFull) || abs(covFull) < eps
        return
    end

    % Subperiod covariance
    n = numel(x);
    bpSub = round(linspace(0, n, nSub + 1));
    covSub = nan(nSub,1);

    for kk = 1:nSub
        idx = bpSub(kk)+1 : bpSub(kk+1);

        xSub = x(idx);
        wSub = w(idx);

        if sum(~(isnan(xSub) | isnan(wSub))) < 3
            covSub(kk) = nan;
        else
            xSub = xSub - mean(xSub, 'omitmissing');
            wSub = wSub - mean(wSub, 'omitmissing');
            covSub(kk) = mean(xSub .* wSub, 'omitmissing');
        end
    end

    covSubMean = mean(covSub, 'omitmissing');

    ssDev = abs((covSubMean - covFull) ./ covFull) .* 100;

end


function itcDev = local_calc_ITCw_deviation(sigma_w, ustar, zeta, z, hc, useCanopyITC)

% w-based integral turbulence characteristics test.
% Uses sigma_w / ustar only.
% Inside canopy, uses Rannik et al.-style canopy parameterization.
% Above canopy, uses Foken-style sigma_w/u* parameterization.

    itcDev = nan;

    if isnan(sigma_w) || isnan(ustar) || ustar <= 0
        return
    end

    measured = sigma_w ./ ustar;

    if isnan(measured) || measured <= 0
        return
    end

    if useCanopyITC && ~isnan(hc) && z <= hc
        model = local_canopy_sigmaw_over_ustar(z, hc);
    else
        model = local_above_canopy_sigmaw_over_ustar(zeta);
    end

    if isnan(model) || model <= 0
        return
    end

    itcDev = abs((model - measured) ./ model) .* 100;

end


function model = local_above_canopy_sigmaw_over_ustar(zeta)

% Foken-style above-canopy sigma_w/u* parameterization.
% From commonly used ITC coefficient table:
%   0 > z/L > -0.032 : c1 = 1.3, c2 = 0
%   -0.032 > z/L     : c1 = 2.0, c2 = 1/8
%
% For stable conditions, a conservative near-neutral fallback is used
% because stable parameterizations are less robust for complex terrain.

    model = nan;

    if isnan(zeta)
        return
    end

    if zeta < -0.032
        model = 2.0 .* abs(zeta).^(1/8);
    elseif zeta < 0
        model = 1.3;
    else
        model = 1.3;
    end

end


function model = local_canopy_sigmaw_over_ustar(z, hc)

% Rannik et al.-style neutral forest-canopy parameterization for sigma_w/u*.
% Coefficients for w:
%   a_i = 1.25
%   alpha_i = 0.9
%   beta_i = 1.2
%   gamma_i = -0.63

    if isnan(z) || isnan(hc) || hc <= 0
        model = nan;
        return
    end

    zh = min(max(z ./ hc, 0), 1);

    ai = 1.25;
    alpha_i = 0.9;
    beta_i = 1.2;
    gamma_i = -0.63;

    model = ai .* ...
        (exp(-alpha_i .* (1 - zh).^beta_i) .* (1 - gamma_i) + gamma_i);

end


function finalFlag = local_dev_to_SSITC_flag(ssDev, itcDev)

% Combined SS + ITC 0-1-2 flag.

    if isnan(ssDev) || isnan(itcDev)
        finalFlag = 2;
    elseif ssDev < 30 && itcDev < 30
        finalFlag = 0;
    elseif ssDev < 100 && itcDev < 100
        finalFlag = 1;
    else
        finalFlag = 2;
    end

end


function finalFlag = local_dev_to_SSonly_flag(ssDev)

% Steady-state-only 0-1-2 flag.
% Intended for slope-flow regimes where ITC assumptions are not valid.

    if isnan(ssDev)
        finalFlag = 2;
    elseif ssDev < 30
        finalFlag = 0;
    elseif ssDev < 100
        finalFlag = 1;
    else
        finalFlag = 2;
    end

end