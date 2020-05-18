 function [out] = runMcStas(p,options,fNames,filename,id)
    
    
    run = tic;
    % Currently this is only supported with MPIRUN due to calls in monitor
    % (WorldSize)
    if options.mpi  > 0 
        McString = ['mpirun'];
        ID = strfind(filename,'.instr');
        filename = [filename(1:ID) 'out'];
        McString = [McString ' ' filename];
        McString = [McString ' --mpi ' num2str(options.mpi) ''];
    else
        McString = ['./'];
        ID = strfind(filename,'.instr');
        filename = [filename(1:ID) 'out'];
        McString = [McString filename];
    end
    
    
    McString = [McString ' --dir ' options.dir '_' num2str(options.generation) '_' num2str(id)];
    McString = [McString ' -n ' num2str(options.ncount)];
    if isfield(options,'seed')
        McString = [McString ' --seed ' num2str(options.seed)];
    end
    McString = [McString ' --gravitation'];
    %McString = [McString ' --overwrite'];
    %McString = [McString ' -c'];
    if ischar(p)
    McString = [McString ' ' p];
    elseif (isstruct(p))
         for i = 1:length(fNames.all)
             McString = [McString ' ' fNames.all{i} '=' eval(['p.' fNames.all{i}]) ];
         end
    end

    try; unix(['rm -rf ' options.dir '_' num2str(options.generation) '_' num2str(id)]);end
 
    [supress,output] = system(McString);
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