function [int_slant, errorFIT, HV] = ...
    cal_virtual_slant(az, el, slant, azV, elV, ZrV, Xp, Yp, levels, minOBSel, HV, steptime, lla)
% cal_virtual_slant estimates virtual slant water vapor from observed GNSS slants.
% It fits an azimuth/elevation model to observed slants, evaluates the model
% along virtual rays, and optionally scales the result with scale-height data.
% az, el    : observed azimuth and elevation in degrees.
% slant     : GNSS slant water vapor values associated with az and el.
% azV, elV  : virtual-ray azimuth and elevation in degrees.
% ZrV       : virtual-ray height where the ray intersects the lateral walls.
% Xp, Yp, levels : 3D tomographic grid.
% minOBSel  : minimum observed elevation used in the fit.
% HV        : scale-height value or directional scale-height array.
% steptime  : assimilation interval used to interpolate HV.
% lla       : longitude, latitude, and orthometric height.

% for testing
% az=az(ind); el=el(ind); slant=slants(ind); PRN=satnum(ind); azV=AZv(indV); elV=ELv(indV);
% ZrV=Cv(indV); steptime=dtimes(t,:); lla=xyzGEO(s,:); HV = HVs(t,:);

mapping = '1/sen';
% For testing, only "1/sen" works when sectors = 1.
% https://benthamopen.com/FULLTEXT/TOASCJ-11-1
%mapping = 'vmf1';         % lla is only used by this mapping.
%mapping = 'chao';
%mapping = 'black&eisner'; % simple geometrical model
%mapping = 'moffett';

top = mean2(levels(:,:,end));
p = 0; % Plot switch.

% Improve the fit by removing observations below "minOBSel".
% Elevations below about 4 degrees usually degrade the result.
idx = find(el < minOBSel);
el(idx) = [];
az(idx) = [];
slant(idx) = [];

% Estimate the model from observed data.
% mfw
if strcmpi(mapping, 'vmf1')
    gpt3GRID = gpt3_1_fast_readGrid();             % Should be read in the main script.
    lon = lla(1); lat = lla(2); alt = lla(3) - 53; % Remove geoid undulation.
    dts = mean(steptime);
    [~,~,~,~,~,~,aw,~,~,~,~,~,~] = gpt3_1_fast(mjuliandate(dts), lat*pi/180, lon*pi/180, alt, 0, gpt3GRID);
    bw = 0.00146;
    cw = 0.04391;
    mfw = (1 + aw .* (1+bw.*(1+cw).^(-1)).^(-1)) ./ (sind(el) + aw.*(sind(el)+bw.*(sind(el)+cw).^(-1)).^(-1));
elseif strcmpi(mapping, 'chao')
    b = 0.00035;
    c = 0.017;
    mfw = 1./(sind(el) + b./(tand(el)+c));
elseif strcmpi(mapping, 'black&eisner')
    a = 1.000009; % Values closer to 1 approach "1/sen".
    mfw = 1./sqrt(1-(cosd(el)./a).^2);
elseif strcmpi(mapping, 'moffett')
    a = 0.05;     % Lower values approach "1/sen".
    mfw = 1./(sind(sqrt(el.^2 + a)));
elseif strcmpi(mapping, '1/sen')
    mfw = 1./sind(el);
end
% a*cosd(az) + b*sind(az) + c/sind(el) + d
const = ones(length(mfw),1);
A = [cosd(az).*mfw, sind(az).*mfw, mfw, const]; % Original model from the paper.
% Y
Y = slant;

% Coefficients from least squares.
coef = (A'*A)\(A'*Y);

% Residuals.
residuos = Y - A*coef;

n = length(Y);
rmse = sqrt(sum(residuos.^2)/n);
mae = 1/n * sum(abs(residuos));
mape = 1/n * sum(abs((residuos)./Y)) * 100;

if strcmpi(mapping, 'vmf1')
    mfw = (1 + aw .* (1+bw.*(1+cw).^(-1)).^(-1)) ./ (sind(elV) + aw.*(sind(elV)+bw.*(sind(elV)+cw).^(-1)).^(-1));
elseif strcmpi(mapping, 'chao')
    mfw = 1./(sind(elV) + b./(tand(elV)+c));
elseif strcmpi(mapping, 'black&eisner')
    mfw = 1./sqrt(1-(cosd(elV)./a).^2);
elseif strcmpi(mapping, 'moffett')
    mfw = 1./(sind(sqrt(elV.^2 + a)));
elseif strcmpi(mapping, '1/sen')
    mfw = 1./sind(elV);
end
int_slant = coef(1).*cosd(azV).*mfw + coef(2).*sind(azV).*mfw + coef(3).*mfw + coef(4);

if p
    xx = 0:360; yy = 0:90;
    [xx, yy] = meshgrid(xx, yy);
    if strcmpi(mapping, 'vmf1')
        mfw = (1 + aw .* (1+bw.*(1+cw).^(-1)).^(-1)) ./ (sind(yy) + aw.*(sind(yy)+bw.*(sind(yy)+cw).^(-1)).^(-1));
    elseif strcmpi(mapping, 'chao')
        mfw = 1./(sind(yy) + b./(tand(yy)+c));
    elseif strcmpi(mapping, 'black&eisner')
        mfw = 1./sqrt(1-(cosd(yy)./a).^2);
    elseif strcmpi(mapping, 'moffett')
        mfw = 1./(sind(sqrt(yy.^2 + a)));
    elseif strcmpi(mapping, '1/sen')
        mfw = 1./sind(yy);
    end
    figure('color', 'w'); hold on
    plane = coef(1).*cosd(xx).*mfw + coef(2).*sind(xx).*mfw + coef(3).*mfw + coef(4);
    surf(xx,yy,plane,'FaceAlpha',0.5); shading interp; c=colorbar;
    scatter3(az, el, slant, 25, slant, 'MarkerEdgeColor','k','MarkerFaceColor','k')
    xlabel('AZ (deg.)'); ylabel('EL (deg.)'); zlabel('ρ_v (kg/m^2)')
    zlim([0, 1100]); xlim([0, 360])
    box on; grid on
    xlabel(c, 'ρ_v (kg/m^2)')
    set(gca, 'fontsize', 12)
    view(-123, 27)
    %export_fig('fit_example.png', '-png', '-m1.5')
    %save('fit_example.mat', 'xx','yy','plane','az','el','slant')
end

% Correction for slants that exit through the domain walls.
% First approach: "A New Unconstrained Approach to GNSS Atmospheric Water Vapor Tomography".
% Zr is the ray height at the domain wall, referenced to the station altitude.
% mul = (1-exp(-(ZrV./HV))) ./ (1-exp(-(top./HV)));
% b = 1;
% a = min(mul);
% xnorm = (b-a).*(mul-a)./(max(mul)-a)+a;
%int_slant = int_slant .* xnorm;
if size(HV, 2) == 1
    int_slant = int_slant .* (1-exp(-(ZrV./HV))) ./ (1-exp(-(top./HV))); % Also compensates for vapor above the domain top.
    %int_slant = int_slant .* (1-exp(-(ZrV.^2./HV.^2))); 
    %int_slant = int_slant .* (1 - (exp(-(ZrV.^2./HV.^2))+exp(-(ZrV./HV)))./2 );     
else
    xSO = Xp(end,1); xNO = Xp(1,1); xSE = Xp(end,end); xNE = Xp(1,end);
    ySO = Yp(end,1); yNO = Yp(1,1); ySE = Yp(end,end); yNE = Yp(1,end);
    xM = mean2(Xp); yM = mean2(Yp);
    % Northern opening angle.
    AB = [xM, yM]-[xNO, yNO];
    CB = [xM, yM]-[xNE, yNE];
    angN = atan2(abs(det([AB;CB])),dot(AB,CB))*180/pi;
    % Eastern opening angle.
    AB = [xM, yM]-[xNE, yNE];
    CB = [xM, yM]-[xSE, ySE];
    angE = atan2(abs(det([AB;CB])),dot(AB,CB))*180/pi;
    % Southern opening angle.
    AB = [xM, yM]-[xSO, ySO];
    CB = [xM, yM]-[xSE, ySE];
    angS = atan2(abs(det([AB;CB])),dot(AB,CB))*180/pi;
    % Western opening angle.
    AB = [xM, yM]-[xNO, yNO];
    CB = [xM, yM]-[xSO, ySO];
    angO = atan2(abs(det([AB;CB])),dot(AB,CB))*180/pi;
    % The sum of the four angles should be 360.
    % angN+angE+angS+angO = 360    
    % N
    indN = azV >= 360-angN/2 | azV < angN/2;
    % E
    indE = azV >= angN/2 & azV < angN/2+angE;    
    % S
    indS = azV >= angN/2+angE & azV < angN/2+angE+angS;       
    % W
    indO = azV >= angN/2+angE+angS & azV < angN/2+angE+angS+angO;       
    % Apply a different HV depending on where virtual rays exit the domain.
    %          1       2     3     4     5
    % ORDER: Center, North, South, East, West.
    int_slant(indN) = int_slant(indN) .* (1-exp(-(ZrV(indN)./HV(2)))) ./ (1-exp(-(top./HV(2)))); 
    int_slant(indE) = int_slant(indE) .* (1-exp(-(ZrV(indE)./HV(4)))) ./ (1-exp(-(top./HV(4))));
    int_slant(indS) = int_slant(indS) .* (1-exp(-(ZrV(indS)./HV(3)))) ./ (1-exp(-(top./HV(3))));
    int_slant(indO) = int_slant(indO) .* (1-exp(-(ZrV(indO)./HV(5)))) ./ (1-exp(-(top./HV(5)))); 
    % Output value.
    HV = mean(HV);
end

errorFIT = [rmse, mae, mape];
return
