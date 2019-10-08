load 'workspace.mat'

%% Determine how many steps is in optimization
if strfind(method,'3step')
    steps=3;
elseif strfind(method,'6step')
    steps=6;
elseif strfind(method,'9step')
    steps=9;
elseif strfind(method,'standard')
    steps=1;
elseif strfind(method,'singleParScan')
    steps=3;
end



if strcmp(printStep,'optimize')
    %% Local file
    %% Open file
    if step == 1
        fileID=fopen('out_CW.txt','w')
        fprintf(fileID,'Optimization Status by CoatingWriter\n')
        fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf(fileID,'Optimizing using the %s method. \n',method)
        fprintf(fileID,'Total steps in process: %i \n',steps)
        fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    else
        % Find last values (crit,I,price)
        fileID=fopen('out_CW.txt','r'); 
        i = 1;
        tline = fgetl(fileID);
        lastFile{i} = tline;
        while ischar(tline)
            i = i+1;
            tline = fgetl(fileID);
            lastFile{i} = tline;
        end
        lastFile{end}=[];
        fclose(fileID); 
        while i > 1
            if strfind(lastFile{i},'Intensity =')
                id = strfind(lastFile{i},'=')
                oldIntensity = eval(lastFile{i}(id+1:end));
                i=1;
            end
            if strfind(lastFile{i},'Price =')
                id = strfind(lastFile{i},'=')
                oldbestPrice = eval(lastFile{i}(id+1:end));
            end
            if strfind(lastFile{i},'Criteria =')
                id = strfind(lastFile{i},'=')
                oldm = eval(lastFile{i}(id+1:end));
            end
            i=i-1;
        end

        fileID=fopen('out_CW.txt','a')
    end

    %% Print details from last iteration
    fprintf(fileID,'\n')
    fprintf(fileID,'\n')
    fprintf(fileID,'Done with step %i of %i \n',step,steps)
    fprintf(fileID,'\t%i variables optimized in %i iterations\n',length(o.parsBest),length(o.parsHistory))

    % Print time
    hours = floor(time / 3600);
    time = time - hours * 3600;
    mins = floor(time / 60);
    secs = time - mins * 60;
    fprintf(fileID,'\tTotal optimization time: %i:%02.0f:%02.0f\n',hours,mins,secs)

    % Print criteria and price
    [m,index]=max(abs(o.criteriaHistory));
    bestPrice=ValueList(index);

    if strcmp(mode,'coating')
        Intensity=o.criteriaBest;
    else
        Intensity=o.criteriaBest*bestPrice;
    end
    fprintf(fileID,'\tValues:\n')
    fprintf(fileID,'\t\tIntensity = %6.4f\n',Intensity)
    fprintf(fileID,'\t\tPrice = %6.2f \n',bestPrice)
    fprintf(fileID,'\t\tCriteria = %d\n',m)
    if step>1
        fprintf(fileID,'\tChange:\n')
        fprintf(fileID,'\t\tIntensity change = %6.2f %%\n',100*(Intensity-oldIntensity)/oldIntensity)
        fprintf(fileID,'\t\tPrice change = %6.2f %% \n',100*(bestPrice-oldbestPrice)/oldbestPrice)
        change=100*(m-oldm)/oldm;	
        if change<0
            sign='-'
        else
            sign='x'
        end
        fprintf(fileID,'\t\tCriteria change = %s%6.2f %%\n',sign,abs(change))
    end
    fclose(fileID)
end


    %% Total status file
    cd ..
    fileID=fopen('RunStatus.txt','w')
    i = 1;
    tline = fgetl(fileID);
    A{i} = tline;
    while ischar(tline)
         i = i+1;
         tline = fgetl(fileID);
         A{i} = tline;
    end
    A(end)=[];
    for i=1:length(A)
        if strcmp(A{i},[instrument_name])
            switch printStep
                case 'optimize' 
                    A{i}=fprintf('%15s |%11s|    %i/%i    |   %5d   | %5d |',instrument_name,'Optimizing',step,steps,Intensity,Price)
                case 'initialize'
                    A{end+1}=fprintf('%15s |%11s|    %i/%i    |   %5s   | %5s |',instrument_name,'Optimizing',0,steps,'n/a','n/a')
                case 'analyze'
                    A{i}=fprintf('%15s |%11s|    %i/%i    |   %5d   | %5d |',instrument_name,'Analysing',step,steps,Intensity,Price)
                case 'finish'
                    A{i}=fprintf('%15s |%11s|    %i/%i    |   %5d   | %5d |',instrument_name,'Finished',step,steps,Intensity,Price)
            end
        done=1;
        end
    end
    for i=1:length(A)
        fprintf(fileID,'%s',A{i})
    end  
fclose(fileID); 
