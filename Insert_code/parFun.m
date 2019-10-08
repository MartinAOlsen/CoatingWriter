function [] = parFun(path,filename)
%PARFUN calls _ifit.m specified in input
%   This function is nessesary for running parallel optimizations, due to
%   the transparacy requirement of the parfor function in matlab parallel
%   computing toolbox

cd(path)
run(filename)
cd ..

end

