 function [out] = runMcStas(p,options,fNames,filename,id)
    
    
    run = tic;
        
    McString = ['mcrun'];
    McString = [McString ' ' filename];
    McString = [McString ' --dir ' options.dir '_' num2str(options.generation) '_' num2str(id)];
    McString = [McString ' --mpi ' num2str(options.mpi) ''];
    McString = [McString ' -n ' num2str(options.ncount)];
    if isfield(options,'seed')
        McString = [McString ' --seed ' num2str(options.seed)];
    end
    McString = [McString ' --gravitation'];
    %McString = [McString ' --overwrite'];
    %McString = [McString ' -c'];
    McString = [McString ' ' p];
    
%     for i = 1:length(fNames.all)
%         McString = [McString ' ' fNames.all{i} '=' eval(['p.' fNames.all{i}]) ];
%     end

    try; unix(['rm -rf ' options.dir '_' num2str(options.generation) '_' num2str(id)]);end
 
    [supress,output] = unix(McString);
    %unix(McString)
    
    
    %% Load data
    if sum(strfind(filename,'analyze')) == 0
        load = LoadCWMonitor([options.dir '_' num2str(options.generation) '_' num2str(id) '/Monitor.txt']);
        out(:,2) = struct2cell(load);
        out(:,1) = fieldnames(load);
        unix(['rm -rf ' options.dir '_' num2str(options.generation) '_' num2str(id)]);
    else
        out=[];
    end
%     out.parameters = p;
    
    
    %fprintf('  Done in %2.2f sec',toc(run))
    %fprintf(' price = %6.0f k€ , intensity = %2.5f  | ',out.price,out.intensity)
    
   
    
end