function [rhov_vox, hgt, voxels, rs_pwv] = rs_voxalization_1d(path2sondes, dt2read, xyCoord, Xp, Yp, levels, verbose)
% rs_voxalization_1d maps radiosonde water-vapor profiles onto tomographic voxels.
% It reads the radiosonde file for one launch time, estimates the sonde path,
% interpolates vapor density to voxel centers, and computes integrated PWV.
% rhov_vox : voxelized water-vapor density.
%      hgt : voxel altitude, approximately the voxel-center mean height.
%   voxels : voxel indices traversed by the sonde [row, column, vertical level].
%   rs_pwv : integrated PWV up to the domain top and total integrated PWV.

% for testing
% path2sondes=p2rs; dt2read=rdtimes(t); xyCoord = [xOBS(t),yOBS(t)]; % In meters.

method = 'linear';
rtime = datevec(dt2read);

% Keep compatibility with older files that do not include minute information.
% Complete timestamp.
l1 = dir([path2sondes,'*_sonde_', num2str(rtime(1)), '-', num2str(rtime(2),'%02i'), '-', num2str(rtime(3),'%02i'), 'T', num2str(rtime(4),'%02i'), num2str(rtime(5),'%02i'), '.dat']);
fn1= [l1.folder,'/',l1.name];
% Timestamp without minutes.
l2 = dir([path2sondes,'*_sonde_', num2str(rtime(1)), '-', num2str(rtime(2),'%02i'), '-', num2str(rtime(3),'%02i'), 'T', num2str(rtime(4),'%02i'), '.dat']);
fn2= [l2.folder,'/',l2.name];

if exist(fn1, 'file') == 2
    af = fn1;
    head = 6;
elseif exist(fn2, 'file') == 2
    af = fn2;
    head = 3;
else
    rhov_vox = nan;
    %rhod_vox = nan;
    voxels = nan;
    hgt = nan;
    rs_pwv = nan;
    %perc = 0;
    if verbose, disp('Radiosonde file does not exist!'); end
    return
end

eps=0.622;
RR=8.31446261815324;
Md=28.9647e-3;
Rd=RR/Md;

text = fileread(af);
newText = splitlines( text );

aux = strsplit(path2sondes, '/');
k = 1;
if ~strcmpi(aux{end-1}, 'smog') % SMOG sondes have a different format.
    PRES = nan(size(newText,1)-head,1);
    HGHT = PRES;
    TEMP = PRES;
    %DWPT = PRES;
    RELH = PRES;
    MIXR = PRES;
    DRCT = PRES;
    SPED = PRES;
    for l = head : size(newText,1)-1
        %-----------------------------------------------------------------------------
        %PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SPED   THTA   THTE   THTV
        %hPa      m      C      C      %   g/kg    deg    m/s      K      K      K
        %-----------------------------------------------------------------------------
        % 1001.4    122   28.6   14.9     43  10.71    170    0.8  301.6  333.5  303.6
        sln = split(newText{l});
        if size(sln,1) > 5 % Complete line.
            PRES(k) = str2num(sln{2})*100;  %#ok<*ST2NM>
            HGHT(k) = str2num(sln{3});
            TEMP(k) = str2num(sln{4})+273.15; % K
            %DWPT(k) = str2num(sln{5});
            RELH(k) = str2num(sln{6})/100;
            MIXR(k) = str2num(sln{7});
            DRCT(k) = str2num(sln{8});
            SPED(k) = str2num(sln{9}); %*0.514444444 % knot to m/s
            k = k + 1;
        end
    end
else
    % SMOG sondes have a different format.
    PRES = nan(size(newText,1)-3,1);
    HGHT = PRES;
    TEMP = PRES;
    %DWPT = PRES;
    RELH = PRES;
    %MIXR = PRES;
    DRCT = PRES;
    SPED = PRES;
    for l = 3 : size(newText,1)-1
        %    Time    Time        Z       p       T      RH    Td      Wspd    Wdir
        %   (min)     (s)      (m)   (hPa)     (C)     (%)    (C)    (m/s)    (deg)
        %       0.      0.    104.  1008.2    23.7     67.    17.2    110.     2.0
        % 1001.4    122   28.6   14.9     43  10.71    170    0.8  301.6  333.5  303.6
        sln = split(newText{l});
        PRES(k) = str2num(sln{5})*100;  %#ok<*ST2NM>
        HGHT(k) = str2num(sln{4});
        TEMP(k) = str2num(sln{6})+273.15; % K
        %DWPT(k) = str2num(sln{5});
        RELH(k) = str2num(sln{7})/100;
        %MIXR(k) = str2num(sln{7});
        DRCT(k) = str2num(sln{9});
        SPED(k) = str2num(sln{8}); %*0.514444444 % knot to m/s
        k = k + 1;
    end
end

idx = find(isnan(HGHT));
PRES(idx) = [];
HGHT(idx) = [];
TEMP(idx) = [];
%DWPT(idx) = [];
RELH(idx) = [];
if ~strcmpi(aux{end-1}, 'smog')
    MIXR(idx) = [];
end
DRCT(idx) = [];
SPED(idx) = [];

% mean ascent rate of 5 m/s
secS = HGHT./5;
uS = SPED.*cosd((270-DRCT));
vS = SPED.*sind((270-DRCT));
xS = zeros(size(secS));
yS = xS;
    
% Sonde launch-site coordinates in meters.
xS(1) = xyCoord(1);
yS(1) = xyCoord(2);

r = 6378137; % m
ntS = length(secS);
for it = 2 : ntS
    xS(it) = xS(it-1)+((0.5*(uS(it)+uS(it-1))*(secS(it)-secS(it-1))))/(r*cosd(yS(1)));
    yS(it) = yS(it-1)+ (0.5*(vS(it)+vS(it-1))*(secS(it)-secS(it-1)))/r;
end
% figure; plot3(xS./1000, yS./1000, HGHT)

if ~strcmpi(aux{end-1}, 'smog')
    rS = MIXR ./ 1000;        % g/kg -> kg/kg
    rhod = PRES ./ (Rd .* TEMP);
    rhov = rhod .* rS;        % kg/m3
else
    % Tetens equation
    % e_sat = 610.78.*exp((17.27.*(TEMP-273.15))./((TEMP-273.15)+273.3)); %P in Pa and T in C
    % The Buck equation
    T1 = 18.678-(TEMP-273.15)./234.5; %P in Pa and T in C
    T2 = (TEMP-273.15)./(257.14+(TEMP-273.15));
    e_sat = 611.21.*exp(T1.*T2);
    eS = RELH.*e_sat;
    rS = eps.*eS./PRES;
    rhod = PRES./(Rd.*TEMP);
    rhov = rhod.*rS;
end

% Total integral.
tpwv = nan(length(HGHT)-1,1);
for k = 1 : length(HGHT)-1
    dhgt = HGHT(k+1) - HGHT(k);
    mrhov= (rhov(k+1) + rhov(k))/2;
    tpwv(k) = dhgt*mrhov;
end
tpwv = sum(tpwv, [], 'omitnan');

% Interpolate to voxel centers.
rhov_vox = nan(size(levels,3)-1, 1);
%rhod_vox = rhov_vox;
hgt = rhov_vox;
voxels = nan(size(levels,3)-1, 3);
kk = 1;
for ii = 1 : size(Xp,2)-1
    for jj = 1 : size(Xp,1)-1
        % old
        %idx_h = xS >= Xp(jj,ii) & xS <= Xp(jj,ii+1) & yS <= Yp(jj,ii) & yS >= Yp(jj+1,ii);
        % new
        % Quadrilateral vertices, clockwise or counterclockwise.
        xv = [Xp(jj,ii), Xp(jj,ii+1), Xp(jj+1,ii+1), Xp(jj+1,ii), Xp(jj,ii)]; % Close the polygon.
        yv = [Yp(jj,ii), Yp(jj,ii+1), Yp(jj+1,ii+1), Yp(jj+1,ii), Yp(jj,ii)];
        % Points to test (xS, yS).
        idx_h = inpolygon(xS, yS, xv, yv);
        
        for k = 1 : size(levels,3)-1
            mli = (levels(jj,ii,k)  +levels(jj+1,ii,k)  +levels(jj+1,ii+1,k)  +levels(jj,ii+1,k))/4;
            mls = (levels(jj,ii,k+1)+levels(jj+1,ii,k+1)+levels(jj+1,ii+1,k+1)+levels(jj,ii+1,k+1))/4;
            ml  = (mli + mls)/2;
            idx_v = HGHT(idx_h) >= mli & HGHT(idx_h) <= mls;
            if ~isempty(idx_v)
                [~,ij] = unique(HGHT);
                rhov_vox(kk) = interp1(HGHT(ij), rhov(ij), ml, method, 'extrap');
                hgt(kk) = ml;
                voxels(kk,1) = jj; % Rows.
                voxels(kk,2) = ii; % Columns.
                voxels(kk,3) = k;  % Levels.
                kk = kk + 1;
            end
        end
    end
end

% If sonde observations end before the domain top.
if HGHT(end) < hgt(end)
    idx = find(hgt > HGHT(end));
    hgt(idx) = nan;
    rhov_vox(idx) = nan;
end

% Integral up to the domain-top altitude using voxelized values.
altMAXpwv = max(max(levels(:,:,end))); % m
tpwv_ = nan(length(hgt)-1,1);
for k = 1 : length(hgt)-1
    dhgt = hgt(k+1) - hgt(k);
    mrhov= (rhov_vox(k+1) + rhov_vox(k))/2;
    tpwv_(k) = dhgt*mrhov;
    if dhgt > altMAXpwv, break; end
end
tpwv_ = sum(tpwv_, [], 'omitnan');

rs_pwv(1) = tpwv_; % Up to the altitude limit using voxelized values.
rs_pwv(2) = tpwv;  % Total using original values.

return
