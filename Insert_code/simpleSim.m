function [monitor] = simpleSim(str,p,options)
    % This is a small function that runs mcstas quickly and only saves the
    % criteria and error. For low ncount this will save a lot of computing
    % time. Designed specifically for the CoatingWriter StepScan mode.


    system(['rm -rf ' options.dir '_stepscan'])

    optionsStr = ['--seed 1 -g -n ' num2str(options.ncount) ' --dir ' options.dir '_stepscan --mpi ' num2str(options.mpi)];
    
    pStr = '';
    F = fieldnames(p)
    for i = 1:length(F)
        pStr = [pStr F{i} '=' eval(['p.' F{i}]) ' '];
    end
    
    [no,output] = system(['mcrun ' str ' ' optionsStr ' ' pStr])
    
    tmp = iLoad([options.dir '_stepscan/Div2d_sample_B.dat'])
    
    monitor.Data.Criteria = sum(sum(tmp.data));
    monitor.Data.Error =  sum(sum(tmp.errors)) / sqrt(size(tmp.errors,1)*size(tmp.errors,2));
end

