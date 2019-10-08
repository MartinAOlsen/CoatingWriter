function mstruct = read_coating(fName)
% reads the output of coatingwriter into a structure

dats = fileread(fName);

% split text into segments

dats = strtrim(strsplit(dats,'\n'));
dats(cellfun(@(C)isempty(C),dats)) = [];
% split on components
type = ismember(cellfun(@(C)C(1:7),dats,'UniformOutput',false),'Segment');

% number of segments
nSeg = str2double(dats{find(type,1,'last')}(8)); %nSeg = sum(type); %

type = cumsum(type);

% stores the segment data
mstruct = cell(1,nSeg);

idx = 1;

for ii = 1:nSeg
    if any(ismember(dats,num2str(ii,'Segment%d')))
        % the segment index ii exists
        mstruct{ii}.str = dats(type==idx);
        % convert string to double
        temp = cell2mat(cellfun(@(C)sscanf(C,'%f,%f,%f ,  %f,  %f,  %f,  %f,'),mstruct{ii}.str(5:end),'UniformOutput',false));
        mstruct{ii}.dat     = temp';
        mstruct{ii}.length  = temp(1,:);
        mstruct{ii}.area    = temp(2,:);
        mstruct{ii}.mLeft   = temp(4,:);
        mstruct{ii}.mRight  = temp(5,:);
        mstruct{ii}.mTop    = temp(6,:);
        mstruct{ii}.mBottom = temp(7,:);

        mstruct{ii}.type = mstruct{ii}.str{3}(19);
        idx = idx+1;
    else
        % no m-value for this segment
        % type = 'N' for these
        mstruct{ii} = struct('str',[],'dat',zeros(0,6),'type','N','length',[]);
    end
end

% resshuffle the the component order to make it compatible with geometry
% data
mstruct = mstruct(end:-1:1);

end