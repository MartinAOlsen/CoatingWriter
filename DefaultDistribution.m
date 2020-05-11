function defaultDistribution=DefaultDistribution(type)

if strcmp(type,'P')
    defaultDistribution='polyN';
    defaultLockSides='all';
elseif strcmp(type,'E')
    defaultDistribution='exp';
    defaultLockSides='all';
elseif strcmp(type,'S')
    defaultDistribution='linear';
    defaultLockSides='all';
elseif strcmp(type,'C')
    defaultDistribution='polyN';
    defaultLockSides='vertical';
elseif strcmp(type,'G')
    defaultDistribution='none';
    defaultLockSides='all';
else
    defaultDistribution='linear';
end
end