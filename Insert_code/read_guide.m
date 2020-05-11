function component = read_guide(fName)
% reads guide geometry from GuideBot/McStas text output into structure
% Made by Sandor, modified by Martin Olsen

dats = fileread(fName);

dats = strsplit(dats,'\n');
dats(cellfun(@(C)isempty(C),dats)) = [];
% split on components
type = ismember(dats,{'P' 'S' 'E' 'G' 'C' 'Sample' 'Moderator'});
type = cumsum(type);

% split up the data
component = cell(1,max(type));

for ii = 1:max(type)
    component{ii}.str  = dats(type==ii);
    component{ii}.type = component{ii}.str{1};
    component{ii}.str  = component{ii}.str(2:end);
    
    % line of sight
    if any(ismember(component{ii}.str,'Los s'))
        component{ii}.los  = 'start';
        component{ii}.str = component{ii}.str(~ismember(component{ii}.str,'Los s'));
    elseif any(ismember(component{ii}.str,'Los e'))
        component{ii}.los  = 'end';
        component{ii}.str = component{ii}.str(~ismember(component{ii}.str,'Los e'));
    else
        component{ii}.los  = 'none';
    end
    
    % read parameters
    for jj = 1:numel(component{ii}.str)
        str0 = component{ii}.str{jj};
        if any(ismember(str0,':'))
            % load the parameters
            str0 = strtrim(strsplit(str0,':'));
            component{ii}.(str0{1}) = str2num(str0{2}); %#ok<ST2NM>
        else
            component{ii}.r = str2num(str0); %#ok<ST2NM>

        end
    end
    
end

end