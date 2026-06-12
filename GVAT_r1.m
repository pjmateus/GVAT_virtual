%
% GVAT - GNSS Water Vapor Tomography r1 (virtual-ray version)
% Main workflow for building the tomographic grid, generating virtual rays,
% assimilating GNSS slant water vapor, and validating against radiosondes.
%
clear; close all; clc

basePATH = 'C:\codex\gvat_github';        % <--- switch the base of a path here
namelistFILE = 'hong_kong_namelist.txt';  % <--- main configuration file
% Enable diagnostic plots. Use with care because many figures can be opened at once.
plt = 0;

disp('--------------------------------------------')
disp('VAT - VApor Tomographic r1 (virtual version)')
disp('--------------------------------------------')
cd(basePATH)

% UNITS: meters and degrees (mm only in siwv reading)
% Control the output of graphs at function level (manual, for testing purposes only)

% Namelist with options.
S = read_namelist(namelistFILE);
names = fieldnames( S );
for i = 1 : length(names)
    eval([names{i},'= S.',names{i},';']);
end

% Fixed directories (must be created by the user)
p2save  = ['solutions/',experiment,'/'];
p2hvs   = ['hvs/',experiment,'/'];  % scale height source
p2exp   = ['experiments/',experiment,'/'];
p2rs    = ['rs/',experiment,'/'];
% Read the HDF5 file.
% Stations outside the tomographic domain should be removed in this section.

if verbose, disp(['Setup for ',upper(experiment),' experiment']); end
% HDF5 path and filename.
d = datevec(d1); yyyy = d(1); mm = d(2);
fn = [p2exp,experiment,'_',num2str(yyyy),'_',num2str(mm,'%02i'),'.h5'];
[sites, xyzGEO] = infoH5file( fn, plt, verbose );

% READ sites data (siwv in kg/m^2)
[idsite, azimuth, elevation, siwv, dt, perc] = readH5file(fn, sites, site2exc, d1, d2, [x0, y0], verbose);

% Some sites may be present in the HDF5 file but have no data in the selected time range.
% Always index xyzCart and sites with the IDs in "est" to select the correct stations.
est = unique(idsite(:,1));
totalNUMsites = length(est);
% Update the station lists.
sites = sites(est,:);
perc  = perc(est);
xyzGEO= xyzGEO(est,:);

% Diagnostic information.
if verbose 
    disp(['nX x nY x nZ = ', num2str(numVx),' x ',num2str(numVy),' x ',num2str(length(resVz))]); disp('Levels:'); disp(resVz);
    idx = elevation == 90;
    disp(['PWV[mean|min|max] = ', num2str(mean(siwv(idx)), '%.2f'),' | ',num2str(min(siwv(idx)), '%.2f'),' | ',num2str(max(siwv(idx)), '%.2f')]); 
end
        
if outliers
    % Check for possible outliers.
    i0 = find(siwv < minIWV); 
    i1 = find(siwv > maxIWV/sind(minEL));
    i2 = find(isnan(siwv));
    i3 = find(isinf(siwv));
    idx = [i0; i1; i2; i3];
    if verbose, if ~isempty(idx), disp(['Removing ',num2str(length(idx)), ' outliers']); end; end
    idsite(idx,:)=[]; azimuth(idx)=[]; elevation(idx)=[]; siwv(idx)=[]; dt(idx)=[];
end

% Create the 3D grid and transform station coordinates with the Gauss transverse Mercator projection.
% NOTES
% 1) The coordinates of the stations must not coincide with the coordinates of the walls
% 2) If the voxels were not rectangles (e.g. voronoi), the formula for calculating the surface has to change
orig = [x0, y0, x0, y0]; % ORDER: see prepare_terrain.m
[Xp, Yp, levels, xyzCart, demSRTM, site_num2del] = ...
    prepare_terrain(orig, rot, numVx, numVy, resVx, resVy, resVz, sites, xyzGEO, est, perc, plt, verbose);
top  = mean2(levels(:,:,end));

% remove stations outside domain
if ~isnan(site_num2del)
    for k = 1 : length(site_num2del)
        disp(['Removing station ',upper(sites(site_num2del(k),:)),'(',num2str(site_num2del(k)),'), it''s outside the domain'])
        idx = find(idsite(:,1) == site_num2del(k));
        idsite(idx,:)=[]; azimuth(idx)=[]; elevation(idx)=[]; siwv(idx)=[]; dt(idx)=[]; 
    end
    % Update the active station list.
    est = unique(idsite(:,1));
    totalNUMsites = length(est);
end

% Enter this block only once after reading the original data.
if ( 1 )    
    % Scale the slants because the station altitude is adjusted to the base of the surface voxel.
    % This has a small impact and compensates for DEM imperfections; it should be improved later.
    oalt = top - xyzGEO(:,3); % The altitude in "xyzGEO" is not updated, only the value in "xyzCart".
    nalt = top - xyzCart(:,3);
    perc = nalt./oalt;
    for k = 1 : totalNUMsites
        ind = idsite(:,1) == est(k);
        siwv(ind) = siwv(ind).*perc(est(k));       
        % Also adjust the station elevation stored in "idsite".
        idsite(ind,4) = xyzCart(est(k),3);
    end
end

% Call once at the beginning to triangulate the domain.
[triD, vertD] = triDomain(Xp, Yp, levels, plt);

% Call once at the beginning to triangulate each voxel.
[tri, vert, indVoxel] = triVoxels(Xp, Yp, levels, plt);
totalNUMvoxels = size(tri,1); % Number of voxels.

% Virtual Rays
% Order: [1:voxelsWALL, 2:fixAZandELcycle, 3:centerOFallVOXELS, 4:riseRAYSinSOMEvoxels, 5:centerOFVOXELSrayMINnum, 6:centerOFVOXELSrayMINnumV2]
% method = [0, 0, 1, 0, 0, 0];
p = [0, 0, 0, 0, 0, 0]; % Enable/disable plots by method for testing only.
[AZv, ELv, virtualID, sitesCOORD, methodID] = virtual_rays(xyzCart, est, Xp, Yp, levels, minEL, method, p, verbose);
% Coordinates of the walls in the tomographic domain of the virtual rays
[Mv, Pv, Cv] = rays(triD, vertD, virtualID, sitesCOORD, AZv, ELv, Xp, Yp, levels, 0, verbose); % these are always the same for any assimilation time

% Lengths of virtual rays inside the voxels for all stations.
% mem
vt.voxelID    = cell(totalNUMsites,1);
vt.distINvoxel= vt.voxelID;
vt.ij         = vt.voxelID; % Contains the indices of removed rays.
vt.indOBS     = vt.voxelID;
vt.wallEXITalt= vt.voxelID;
vt.numINI     = vt.voxelID;
vt.numDEL     = vt.voxelID;
vt.numEND     = vt.voxelID;
totalNUMvirtualRays = 0;
tStart = tic;
for s = 1 : totalNUMsites
    % Select all elements for the same station.
    ind = virtualID == est(s);
    xyzCi = xyzCart(est(s),:);
    xyzCv = [Mv(ind), Pv(ind), Cv(ind)];
    tnv = sum(ind);
    voxelIDvt = nan(totalNUMvoxels, tnv); DISTvt = nan(totalNUMvoxels, tnv); diffvt = nan(tnv, 1); wallEXITalt = diffvt;
    elev = ELv(ind);
    azim = AZv(ind);
    mid  = methodID(ind);
    parfor k = 1 : tnv
        [voxelIDvt(:,k), DISTvt(:,k), diffvt(k)] = inVoxels(tri, vert, indVoxel, xyzCi, xyzCv(k,:));         
        % Height difference between the station and the ray when it reaches the wall.
        % "wallRAYSexitAJUST" compensates for the 10 m extension added in triDomain.m.
        distINC = norm(xyzCv(k,:) - xyzCi) - wallRAYSexitAJUST;
        wallEXITalt(k) = sind(elev(k)) * distINC; 
    end
    i0 = find( diffvt > 1 | isnan(diffvt) | isinf(diffvt) ); % Remove rays that pass below the surface, e.g., site ROCA vs. Serra de Sintra.
    sDISTvt = sum(DISTvt, 'omitnan');
    i1 = find( sDISTvt < delRAYSlessDIST )'; % Remove rays shorter than "delRAYSlessDIST" in meters.
    i2 = find( sDISTvt > delRAYSmoreDIST )'; % Remove rays longer than "delRAYSmoreDIST" in meters.
    ij = unique([i0; i1; i2]);
    voxelIDvt(:,ij)= [];
    DISTvt(:,ij)   = [];
    wallEXITalt(ij)= [];
    indOBS = find(virtualID == est(s)); indOBS(ij) = [];
    % Avoid storing so many NaN values.
    tdel = find( all(isnan(voxelIDvt),2) );
    voxelIDvt(tdel,:)= [];
    DISTvt(tdel,:)   = [];
    % Save in the structure for use in the main loop.
    vt.voxelID{s}    = voxelIDvt;
    vt.distINvoxel{s}= DISTvt;
    vt.ij{s}         = ij;
    vt.wallEXITalt{s}= wallEXITalt;
    vt.indOBS{s}     = indOBS; 
    vt.numINI{s}     = size(diffvt,1);
    vt.numDEL{s}     = length(ij);
    vt.numEND{s}     = vt.numINI{s} - vt.numDEL{s};
    totalNUMvirtualRays = totalNUMvirtualRays + vt.numEND{s};
    % info
    if verbose
        fprintf('site: %s(%02i) | virtual rays[ini|end|del]=%5i |%5i |%4i | rays[< %im, del]= %3i | rays[> %im, del]= %3i\n',upper(sites(est(s),:)),est(s),vt.numINI{s},vt.numEND{s},vt.numDEL{s},delRAYSlessDIST,length(i1),delRAYSmoreDIST,length(i2))
    end
end
tEnd = toc(tStart);
disp(['... done in ',num2str(tEnd, '%.2f'),' sec '])

% Radiosonde data.
if runRStimes
    % vector with radiosonde instants (faster, only for reading dates)
    [dtimes, rdtimes, lonRS, latRS] = read_rs_datetimes(p2rs, assimilationTIME);
    idx = find(dtimes(:,2) < d1 | dtimes(:,1) > d2);
    dtimes(idx,:) = [];
    rdtimes(idx,:) = []; % Correct sonde timing so validation does not need to call "dir" repeatedly.
    lonRS(idx) = []; latRS(idx) = [];
    xOBS = nan(size(lonRS,1),1); yOBS = xOBS;
    for k = 1 : size(lonRS,1)
        [xOBS(k),yOBS(k)] = gauss(x0, y0, lonRS(k), latRS(k));
    end
else
    error(['GVAT r1 only works for the RS times (located in ',p2rs,')'])
end
totalTIMES = size(dtimes,1);

% Total number of different designs
activeSITES = zeros(totalTIMES, totalNUMsites);
indDT = cell(totalTIMES, 1); 
for t = 1 : totalTIMES
    % Assimilation range.
    i0 = dtimes(t,1); i1 = dtimes(t,2);
    ind0 = find(dt >= i0 & dt <= i1);
    if ~isempty(ind0)
        el_ = elevation(ind0);
        for s = 1 : totalNUMsites
            ind1 = find(idsite(ind0,1) == est(s));            
            ind2 = find(el_(ind1) ~= 90); % Slants are mandatory.
            ind3 = find(el_(ind1) == 90); % Vertical observations are mandatory.
            if ~isempty(ind1) && ~isempty(ind2) && ~isempty(ind3) 
                activeSITES(t,s) = 1;                
            end
        end
    end
    indDT{t} = ind0; % Used in the main loop.
end
[uniqueDESIGNS,~,idDESIGNS] = unique(activeSITES, 'rows', 'stable');
if verbose
    disp([num2str(size(uniqueDESIGNS,1)),' distinct designs found']); 
    disp(uniqueDESIGNS)
end
% Making unique designs
tStart = tic;
if verbose, disp('Making unique designs'); end
numDESIGNS = size(uniqueDESIGNS,1); 
D = cell(numDESIGNS,1);
for k = 1 : numDESIGNS
    k1 = 0;
    numLL = 0;
    for s = 1 : totalNUMsites
        if uniqueDESIGNS(k,s) == 1, numLL = numLL + vt.numEND{s}; end 
    end
    dg = zeros(numLL, totalNUMvoxels);
    for s = 1 : totalNUMsites
        if uniqueDESIGNS(k,s) == 1
            for k0 = 1 : vt.numEND{s}
                id = vt.voxelID{s}(:,k0);
                d  = vt.distINvoxel{s}(:,k0);
                dg(k1+k0, id(~isnan(id))) = d(~isnan(d));
            end
            k1 = k1 + vt.numEND{s};
        end
    end
    D{k} = dg;
end
tEnd = toc(tStart);
if verbose, disp(['... done in ',num2str(tEnd, '%.2f'),' sec ']); end
clearvars dg numLL uniqueDESIGNS el_ azim elev
            
% Scale height.
[Hvdt, HVs] = getScaleHeight(experiment, p2hvs, hv_source, hv_value, dtimes, verbose);
n = size(Xp,1)-1;
m = size(Xp,2)-1;
l = size(levels,3)-1;

% FIXdata
FIXdata.paths = cell(3,1); 
FIXdata.paths{1}=p2save; FIXdata.paths{2}=p2hvs; FIXdata.paths{3}=p2rs;
FIXdata.assimilationTIME = assimilationTIME;
FIXdata.Xp = Xp; FIXdata.Yp = Yp; FIXdata.levels = levels;
% Coordinates of voxel centers.
Xm = nan(n,m,l); Ym = Xm; Zm = Xm;
for kk = 1 : l
    for jj = 1 : m
        for ii = 1 : n
            Xm(ii,jj,kk) = mean([Xp(ii,jj), Xp(ii+1,jj), Xp(ii+1,jj+1), Xp(ii,jj+1)]);
            Ym(ii,jj,kk) = mean([Yp(ii,jj), Yp(ii+1,jj), Yp(ii+1,jj+1), Yp(ii,jj+1)]);
            Zm(ii,jj,kk) = mean([levels(ii,jj,kk),  levels(ii+1,jj,kk),  levels(ii+1,jj+1,kk),  levels(ii,jj+1,kk), levels(ii,jj,kk+1),levels(ii+1,jj,kk+1),levels(ii+1,jj+1,kk+1),levels(ii,jj+1,kk+1)]);
        end
    end
end
FIXdata.Xm = Xm; FIXdata.Ym = Ym; FIXdata.Zm = Zm;
FIXdata.times          = [dtimes(:,1), (dtimes(:,1)+dtimes(:,2))./2, dtimes(:,2)]; % Time intervals.
FIXdata.rtimes         = rdtimes;  % Correct observation times.
FIXdata.coordSITEScart = xyzCart(est,:);
FIXdata.coordSITESgeo  = xyzGEO(est,:); % Keeps the original station height and removes stations outside the domain.
FIXdata.nameSITES      = sites(est,:);
FIXdata.ll0            = [x0, y0]; % System origin in degrees for converting WRF and RS data during validation.
% Save station PWV for validation.
for k = 1 : totalNUMsites
    ind = (elevation == 90) & idsite(:,1) == est(k);
    FIXdata.obsIWV{k}   = siwv(ind);
    FIXdata.obsIWVdt{k} = dt(ind);
end

% TOMOdata
TOMOdata.x = cell(totalTIMES,1);
TOMOdata.raysINvoxels= TOMOdata.x;
TOMOdata.totalSlants = zeros(totalTIMES,1);
TOMOdata.numOFsites  = nan(totalTIMES,1);
TOMOdata.cond = TOMOdata.numOFsites;
TOMOdata.rank = TOMOdata.numOFsites;
TOMOdata.fit  = nan(totalTIMES,3); 
TOMOdata.minmaxELobs = nan(totalTIMES,3);
TOMOdata.rmse = nan(totalTIMES,1);
TOMOdata.Hv   = nan(totalTIMES,1);

disp('*** Only virtual rays will be added to the design ***');
tStart0 = tic;
mgMISSINGsites = cell(totalTIMES,1);
d0m = zeros(totalTIMES,1); d1m = d0m;
rmseRS = nan(totalTIMES,1); 
biasRS = rmseRS; 
corrRS = rmseRS; 

% Instant loop.
for t = 1 : totalTIMES
    
    % Assimilation range.
    i0 = dtimes(t,1);
    i1 = dtimes(t,2);
    
    % All observations at instant "t".
    obsIND= indDT{t};    
    az    = azimuth(obsIND); 
    el    = elevation(obsIND);
    slants= siwv(obsIND);
    
    % Check whether GNSS data are available at instant "t".
    if sum(activeSITES(t, :)) > 0
        
        % Counters, memory allocation, and diagnostics.
        % Instant-level diagnostics.
        if verbose && atINSTANTlevel
            tStart1 = tic; 
            idx = logical(activeSITES(t, :)); 
        end
        k0 = 1;
        errorFITall = nan(totalNUMsites,3);        
        missingSite = 'site(s):';        
        A = D{idDESIGNS(t)};
        b = nan(size(A, 1), 1);        
        
        % Site loop.
        for s = 1 : totalNUMsites
            % Check whether site "s" has data at instant "t".
            if activeSITES(t, s) == 1
                
                % Site-level diagnostics.
                if verbose && atSITElevel, tStart2 = tic; end
                % Observation indices for site "s".
                ind = find(idsite(obsIND,1) == est(s)); 
                %figure(888); r = 90 - el(ind); polarplot(deg2rad(az(ind)), r, 'o');
                % Virtual-ray indices for site "s".
                indV = vt.indOBS{s};               

                % virtual slants
                [slantvt, errorfit, Hv] = ...
                    cal_virtual_slant(az(ind), el(ind), slants(ind), AZv(indV), ELv(indV), vt.wallEXITalt{s}, Xp, Yp, levels, minOBSel, HVs(t,:), dtimes(t,:), xyzGEO(est(s),:));
                errorFITall(s,:) = mean(errorfit, 1, 'omitnan'); 
                      
                % Fill array b.
                k1 = k0 + vt.numEND{s} - 1;
                b(k0:k1) = slantvt;
                k0 = k0 + vt.numEND{s};
                % Site-level verbose output.
                if verbose && atSITElevel
                    tEnd2 = toc(tStart2);
                    fprintf('site: %s(%02i) | rays[obs|vir|del]=%5i |%5i |%4i | Hv= %6.1f | fit rmse=%4.1f | step=%4i(%4i) | t= %4.2fs\n',...
                        upper(sites(est(s),:)),est(s),length(ind),vt.numEND{s},vt.numDEL{s},Hv,errorFITall(s,1),t,totalTIMES,tEnd2)
                end
            else
                % Add information about missing sites at instant "t".
                missingSite = [missingSite,' ',upper(sites(est(s),:))];  %#ok<AGROW>
            end
        end
        mgMISSINGsites{t} = missingSite; d0m(t) = i0; d1m(t) = i1;                      

        % System inversion using NNLS or lsqminnorm.
        [TOMOdata.x{t}, TOMOdata.rmse(t,1)] = ...
            g_inverse(A, b, avoid_neg, accy, quadratic_optimization);

        % Validation is time-consuming but can be reused by plot_validation_data.m.
        % RMSE against models is the 3D RMSE.
        [rmseRS(t),biasRS(t),corrRS(t),TOMOdata.tomo_rs_rhov_1d{t},TOMOdata.rs_rhov_1d{t},TOMOdata.rs_hgt_1d{t}] = ...
            validation_with_rs(Xp, Yp, levels, rdtimes(t), TOMOdata.x{t}, p2rs, [xOBS(t),yOBS(t)], plt);

        TOMOdata.cond(t) = cond(A); 
        TOMOdata.rank(t) = rank(A);
        % to determine the number of rays passing through each voxel, note that we are changing "A"
        A(A > 0) = 1;
        TOMOdata.raysINvoxels{t} = sum(A,1);
        TOMOdata.numOFsites(t) = sum(activeSITES(t, :));
        TOMOdata.totalSlants(t) = length(obsIND);
        TOMOdata.fit(t,:) = mean(errorFITall, 1, 'omitnan');
        % info about the lowest, fixed and highest elevation observed at time "t"
        TOMOdata.minmaxELobs(t,1) = min(el);
        TOMOdata.minmaxELobs(t,2) = minOBSel; % fix
        TOMOdata.minmaxELobs(t,3) = max(el(el~=90));
        TOMOdata.Hv(t) = Hv;
        
        % Instant-level verbose output.
        if verbose && atINSTANTlevel
            tEnd1 = toc(tStart1);                
            vald = ['rmse[rs]= ',num2str(rmseRS(t),'%.1f')];
            fprintf('step=%4i(%4i) | rays[obs|vir|del]=%6i |%6i |%4i | Hv= %6.1f | fit rmse=%4.1f | lsm rmse=%4.1f | %s | D=%2i(%2i) | t= %4.2fs\n', ...
                t,totalTIMES,TOMOdata.totalSlants(t),sum([vt.numEND{idx}]),sum([vt.numDEL{idx}]),Hv,mean(errorFITall(:,1),'omitnan'),TOMOdata.rmse(t,1),vald,idDESIGNS(t),max(idDESIGNS),tEnd1);
        end
    else
        if verbose, disp(['No data between ', datestr(i0,'dd/mm/yy HH:MM:SS'), ' and ',datestr(i1,'dd/mm/yy HH:MM:SS')]); end
    end  
end
for t = 1 : totalTIMES
    if length(mgMISSINGsites{t}) >= 9
        disp([mgMISSINGsites{t}, ' | no data between ', datestr(d0m(t),'dd/mm/yy HH:MM:SS'), ' and ',datestr(d1m(t),'dd/mm/yy HH:MM:SS'), ' | step: ',num2str(t)]);
    end
end
tEnd0 = toc(tStart0);                
vald = ['mean rmse[rs]= ',num2str(mean(rmseRS,'omitnan'),'%.2f')];
fprintf('mean fit rmse= %.1f | mean lsm rmse= %.1f | %s | mean slants= %i | mean Cond.= %.6f\n',...
    mean(TOMOdata.fit(:,1),'omitnan'), mean(TOMOdata.rmse,'omitnan'), vald, mean(TOMOdata.totalSlants,'omitnan'), mean(TOMOdata.cond,'omitnan'))
disp(['... all done in ',num2str(tEnd0/60, '%.2f'),' min '])

% ************************************ END CORE

if saveSOLUTION
    % validation data
    rmseRS(idxnan) = []; biasRS(idxnan) = []; corrRS(idxnan) = [];
    TOMOdata.rmseRS  = rmseRS;  TOMOdata.biasRS  = biasRS;  TOMOdata.corrRS  = corrRS;
    
    FNversion = ['_minEL',num2str(minEL,'%02i'),'_a',num2str(assimilationTIME,'%02i'),'min_x',num2str(numVx,'%02i'),'y',num2str(numVy,'%02i'),...
    'z',num2str(length(resVz),'%02i'),'_rx',fix(num2str(resVx/1000)),'km_ry',fix(num2str(resVy/1000)),'km_hv',...
    upper(replace(hv_source,'/','_')),'_rot',num2str(rot),'_v',num2str(version)];

    disp(['Saving version: ',FNversion])
    save([p2save,experiment,FNversion,'.mat'], 'TOMOdata', 'FIXdata', 'S', '-v7.3')
end
