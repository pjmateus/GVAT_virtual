function [Hvdt, HV] = getScaleHeight(experiment, p2hvs, hv_source, hv_value, dtimes, verbose)
% getScaleHeight returns scale-height values for each assimilation interval.
% Values are either constant or linearly interpolated from an external file,
% avoiding repeated interpolation inside the main tomographic loop.
% Scale height is the altitude increase required for a quantity to decrease by a factor of 1/e.

if strcmpi(hv_source, 'const')
    % In this case, "ajustHV" is not applied.
    Hvdt = mean(dtimes,2);
    HV = zeros(size(dtimes,1),1) + hv_value;
    return
else
    if verbose, disp(['Reading scale height values from hv_',hv_source,'_',experiment,'.txt file']); end
    T = table2array(readtable([p2hvs,'hv_',hv_source,'_',experiment,'.txt']));
    if size(T,2) == 5      % Hourly information.
        Hvdt = datenum(T(:,1),T(:,2),T(:,3),T(:,4),zeros(size(T,1),1),zeros(size(T,1),1));
        Hv   = T(:,5);
    elseif size(T,2) == 6  % Minute-level information.
        Hvdt = datenum(T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),zeros(size(T,1),1));
        Hv   = T(:,6);       
    end   
    if verbose, disp(['First Hv date in file: ',datestr(Hvdt(1),'dd/mm/yy HH:MM:SS'),' | last date: ',datestr(Hvdt(end),'dd/mm/yy HH:MM:SS')]); end

    % Linear interpolation for the "dtimes" instants.
    dts = mean(dtimes,2);
    idx = dts < Hvdt(1) | dts > Hvdt(end);
    if sum(idx) ~= 0
        disp('Warning: extrapolation of scale height (Hv) !')
    end

    if size(Hv,2) == 1
        HV = interp1(Hvdt, Hv, dts, 'linear', 'extrap');
    elseif size(Hv,2) == 5 % Hv order: Center, North, South, East, West.
        HV = nan(length(dts),5);
        for k = 1 : 5
            HV(:,k) = interp1(Hvdt, Hv(:,k), dts, 'linear', 'extrap');
        end
    else
        error('Hv must have size n-by-1 or n-by-5')
    end

    idx = find(HV < 100);
    if ~isempty(idx)
        disp('Hv values below 100 were found; using the default value of 100')
        HV(idx) = 100;
    end
end

return
