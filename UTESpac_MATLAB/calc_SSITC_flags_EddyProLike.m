function [fluxQC_row, fluxQCHeader] = calc_SSITC_flags_EddyProLike( ...
    fluxQC_row, fluxQCHeader, ...
    jj, ii, sonHeight, info, ...
    wPF_P, uPF_P, vPF_P, Theta_v_sonP, ...
    H2Op, rhov_externalp, rho_CO2p, rhoc_externalp, ...
    tau, L, ...
    rotatedSonFlag, TsonFlag, H2OFlag, CO2Flag)

% calc_SSITC_flags_EddyProLike
%
% Computes EddyPro-like Mauder and Foken 2004 / CarboEurope-style
% 0/1/2 SSITC flags:
%
%   0 = high quality
%   1 = moderate quality / usable
%   2 = poor quality / discard or gap-fill
%
% Steady-state test:
%   compares full-period covariance to mean of 5-min subperiod covariances.
%
% ITC test:
%   uses sigma_w/u* only, following the simplified CarboEurope-style
%   0/1/2 SSITC framework.

    numFluxQCVariables = 4;

    % -----------------------------
    % Settings
    % -----------------------------
    SSITC_subAvgMin = 5;

    if isfield(info,'SSITC_subAvgMin')
        SSITC_subAvgMin = info.SSITC_subAvgMin;
    end

    if mod(info.avgPer, SSITC_subAvgMin) ~= 0
        error('info.avgPer must be divisible by info.SSITC_subAvgMin.');
    end

    nSub = info.avgPer / SSITC_subAvgMin;

    if isfield(info,'displacementHeight')
        d = info.displacementHeight;
    else
        d = 0;
    end

    z_eff = max(sonHeight - d, 0.1);

    % -----------------------------
    % ITC_w calculation
    % -----------------------------
    numTauVariables = 14;
    tau_col = 3 + (ii-1)*numTauVariables;

    if size(tau,2) >= tau_col && ~isnan(tau(jj,tau_col)) && tau(jj,tau_col) >= 0
        ustar = sqrt(tau(jj,tau_col));
    else
        ustar = nan;
    end

    sigma_w = std(wPF_P, 'omitmissing');

    L_col = 2 + (ii-1);

    if size(L,2) >= L_col && ~isnan(L(jj,L_col)) && L(jj,L_col) ~= 0
        zeta = z_eff ./ L(jj,L_col);
    else
        zeta = nan;
    end

    itcWDev = local_calc_ITCw_deviation(sigma_w, ustar, zeta);

    % -----------------------------
    % TAU_SSITC_TEST
    % -----------------------------
    tauSSDev_uw = local_calc_SS_deviation(uPF_P, wPF_P, nSub);
    tauSSDev_vw = local_calc_SS_deviation(vPF_P, wPF_P, nSub);

    tauSSDev = max([tauSSDev_uw, tauSSDev_vw], [], 'omitmissing');

    tauFlag = local_combine_SSITC_flag(tauSSDev, itcWDev);

    if rotatedSonFlag(jj)
        tauFlag = 2;
    end

    fluxQC_row(2+(ii-1)*numFluxQCVariables) = tauFlag;
    fluxQCHeader{2+(ii-1)*numFluxQCVariables} = ...
        strcat(num2str(sonHeight),'m:TAU_SSITC_TEST');

    % -----------------------------
    % H_SSITC_TEST
    % -----------------------------
    HSSDev = local_calc_SS_deviation(wPF_P, Theta_v_sonP, nSub);

    HFlag = local_combine_SSITC_flag(HSSDev, itcWDev);

    if rotatedSonFlag(jj) || TsonFlag(jj)
        HFlag = 2;
    end

    fluxQC_row(3+(ii-1)*numFluxQCVariables) = HFlag;
    fluxQCHeader{3+(ii-1)*numFluxQCVariables} = ...
        strcat(num2str(sonHeight),'m:H_SSITC_TEST');

    % -----------------------------
    % LE_SSITC_TEST
    % -----------------------------
    if ~isempty(H2Op) && ~isempty(rhov_externalp)
        H2O_WPL_P = H2Op + rhov_externalp .* 1e3;  % g/m^3

        LESSDev = local_calc_SS_deviation(wPF_P, H2O_WPL_P, nSub);
        LEFlag = local_combine_SSITC_flag(LESSDev, itcWDev);

        if rotatedSonFlag(jj) || H2OFlag(jj)
            LEFlag = 2;
        end
    else
        LEFlag = nan;
    end

    fluxQC_row(4+(ii-1)*numFluxQCVariables) = LEFlag;
    fluxQCHeader{4+(ii-1)*numFluxQCVariables} = ...
        strcat(num2str(sonHeight),'m:LE_SSITC_TEST');

    % -----------------------------
    % FC_SSITC_TEST
    % -----------------------------
    if ~isempty(rho_CO2p) && ~isempty(rhoc_externalp)
        CO2_WPL_P = rho_CO2p + rhoc_externalp;  % kg/m^3

        FCSSDev = local_calc_SS_deviation(wPF_P, CO2_WPL_P, nSub);
        FCFlag = local_combine_SSITC_flag(FCSSDev, itcWDev);

        if rotatedSonFlag(jj) || CO2Flag(jj)
            FCFlag = 2;
        end
    else
        FCFlag = nan;
    end

    fluxQC_row(5+(ii-1)*numFluxQCVariables) = FCFlag;
    fluxQCHeader{5+(ii-1)*numFluxQCVariables} = ...
        strcat(num2str(sonHeight),'m:FC_SSITC_TEST');

end


function ssDev = local_calc_SS_deviation(xPrime, yPrime, nSub)

    ssDev = nan;

    if isempty(xPrime) || isempty(yPrime)
        return
    end

    if numel(xPrime) ~= numel(yPrime)
        return
    end

    valid = ~(isnan(xPrime) | isnan(yPrime));

    if sum(valid) < 10
        return
    end

    covFull = mean(xPrime .* yPrime, 'omitmissing');

    if isnan(covFull) || abs(covFull) < eps
        return
    end

    n = numel(xPrime);
    bpSub = round(linspace(0, n, nSub + 1));
    covSub = nan(nSub,1);

    for kk = 1:nSub
        idx = bpSub(kk)+1 : bpSub(kk+1);
        covSub(kk) = mean(xPrime(idx) .* yPrime(idx), 'omitmissing');
    end

    covSubMean = mean(covSub, 'omitmissing');

    ssDev = abs((covSubMean - covFull) ./ covFull) .* 100;

end


function itcDev = local_calc_ITCw_deviation(sigma_w, ustar, zeta)

    itcDev = nan;

    if isnan(sigma_w) || isnan(ustar) || isnan(zeta) || ustar <= 0
        return
    end

    itcMeasured = sigma_w ./ ustar;

    % EddyPro-like / Foken-style simplified reference behavior.
    % This is not guaranteed to be bitwise identical to EddyPro, but follows
    % the common sigma_w/u* ITC logic.
    if zeta < 0
        itcModel = 1.25 .* (1 - 3 .* zeta).^(1/3);
    else
        itcModel = 1.25;
    end

    if isnan(itcModel) || itcModel <= 0
        return
    end

    itcDev = abs((itcMeasured - itcModel) ./ itcModel) .* 100;

end


function qFlag = local_combine_SSITC_flag(ssDev, itcDev)

    if isnan(ssDev) || isnan(itcDev)
        qFlag = 2;
    elseif ssDev < 30 && itcDev < 30
        qFlag = 0;
    elseif ssDev < 100 && itcDev < 100
        qFlag = 1;
    else
        qFlag = 2;
    end

end