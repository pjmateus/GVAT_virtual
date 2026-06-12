function S = read_namelist( filename )
% read_namelist parses a WRF-style namelist into a MATLAB structure.
% Each namelist variable is returned as a field in S.
% S = READ_NAMELIST(filename) returns values from WRF-type namelist
%   This function reads input namelist and returns the values as
%   the fields of a structure. 
%

% Open a text file
if ischar(filename)
    try
        fid=fopen(filename,'r');
    catch
        disp(['ERROR: namelist file ' filename ' not found!']);
        S = nan;
        return
    end
else
    fid=filename;
end

while ~feof(fid)
    line = fgetl(fid);
    if ~strncmpi(strtrim(line),'!',1) && length(line) > 1
        if strcmp(line(end), ','), line(end)=''; end
        % The line does not start with a comment.
        sline = split(line);
        var = strtrim(sline{1});
        value = strtrim(sline{3});
        i = 4;
        while endsWith(value,',')
            value = [value, strtrim(sline{i})];
            i = i + 1;
        end
        if strcmp(value(end), ';'), value(end)=''; end
        evalc(['S.',var,'=',value,';']);
    end
end
fclose(fid);

return
