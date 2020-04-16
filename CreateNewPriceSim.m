function  []=CreateNewPriceSim(sourceDir,priceList)

%% Initial variables
ParWorks=1;
    % path split: 
[upperPath, deepestFolder, ~] = fileparts(sourceDir);

    % cd out from folder


    % determine if OS uses / or \
nslash=length(strfind(upperPath,'/'));
nbackslash=length(strfind(upperPath,'\'));
if nslash<nbackslash
    slash='\';
else
    slash='/';
end

    % launch files first lines:
    % Cluster
compile{1}=[];
compilePY{1}=[];
launch_cluster{1}='#!/bin/bash';
launch_cluster{end+1}='# Made by CoatingWriter for guide_bot';
launch_cluster{end+1}='# Replaces a former guide_bot file';
launch_cluster{end+1}='# Will launch all the batch files (works on DMSC cluster)';
launch_cluster{end+1}='';
launch_cluster{end+1}='thisFolder = cd;';
launch_cluster{end+1}='';
launch_cluster{end+1}='cd brilliance_refference';
launch_cluster{end+1}='sbatch brilliance.batch';
launch_cluster{end+1}='cd(thisFolder)'


RunAll{1}='#!/bin/bash';
RunAll{end+1}='# Made by CoatingWriter for guide_bot';
RunAll{end+1}='# Replaces a former guide_bot file';
RunAll{end+1}='# Will launch all the batch files (works on DMSC cluster)';
RunAll{end+1}='';
RunAll{end+1}='thisFolder = cd;';
RunAll{end+1}='';
RunAll{end+1}='cd brilliance_refference';
RunAll{end+1}='mcrun-py -c -n 0 --mpi=1 brilliance.instr ';
RunAll{end+1}='sbatch brilliance.batch';
RunAll{end+1}='cd(thisFolder)'

    % Home
launch_home{1}='% This script will run all optimizations in price-scan one at a time';
if ParWorks==1
    launch_home{end+1}='% A significant performance increase is possible if "parallel computing toolbox" is installed';
    launch_home{end+1}='% To run in parallel, use the launch_home_parallel.m file instead';
end
launch_home{end+1}='';
launch_home{end+1}='%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
launch_home{end+1}='% This script is generated using CoatingWriter - contact Martin A. Olsen (martinolsen.mo@gmail.com) for help';
launch_home{end+1}='%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
launch_home{end+1}='';
launch_home{end+1}='thisFolder = cd;';
launch_home{end+1}='';
launch_home{end+1}='cd brilliance_refference';
launch_home{end+1}='run(''brilliance_ifit.m'');';
launch_home{end+1}='cd(thisFolder)';

    % Home with parallel computing toolbox
launch_home_parallel{1}='% This script will run all optimizations in price-scan in parallel';
launch_home_parallel{end+1}='% If matlab parallel computing toolbox is not installed, this will fail';
launch_home_parallel{end+1}='% If problems occour while running in parallel, try using the non-parallel version, "launch_home.m"';
launch_home_parallel{end+1}='';
launch_home_parallel{end+1}='%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
launch_home_parallel{end+1}='% This script is generated using CoatingWriter - contact Martin A. Olsen (martinolsen.mo@gmail.com) for help';
launch_home_parallel{end+1}='%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
launch_home_parallel{end+1}='';
launch_home_parallel{end+1}='Folder=cd;';
launch_home_parallel{end+1}='[upperPath, deepestFolder, ~] = fileparts(Folder);';
launch_home_parallel{end+1}='% determine if OS uses / or \';
launch_home_parallel{end+1}=['nslash=length(strfind(upperPath,''/''));'];
launch_home_parallel{end+1}=['nbackslash=length(strfind(upperPath,''\''));'];
launch_home_parallel{end+1}='if nslash<nbackslash';
launch_home_parallel{end+1}=['    slash=''\'';'];
launch_home_parallel{end+1}='else';
launch_home_parallel{end+1}=['    slash=''/'';'];
launch_home_parallel{end+1}='end';
launch_home_parallel{end+1}='';
launch_home_parallel{end+1}='';
launch_home_parallel{end+1}='% List of all paths and files:'
launch_home_parallel{end+1}=''   
launch_home_parallel{end+1}=['fileList{1}=''brilliance_ifit.m'';']       
launch_home_parallel{end+1}=['pathList{1}=[Folder slash ''brilliance_refference''];']                



for i=1:length(priceList)
    fprintf('Creating folder for optimization with price: %2.2f\n',priceList(i))
    thisPrice=priceList(i);
    %% Dublicate folder and change price
    copyfile([upperPath slash deepestFolder],[upperPath slash deepestFolder '_price_' num2str(thisPrice)]);
    cd([upperPath slash deepestFolder '_price_' num2str(thisPrice)]);
    
    % List all _ifit.m files in dir:
    s = dir('*_ifit.m');
    file_list = {s.name};
    ifits={};
    
    for j=1:length(file_list)
        tmp=char(file_list(j));
        if strcmp(tmp(1:3),'run')==0
            if numel(ifits)==0
                ifits{1}={tmp};
            else
                ifits{end+1}={tmp};
            end
        else
            {};
        end
    end
    instrumentName=char(ifits{1});
    instrumentName = strtok(instrumentName, '_');
    
    % Give warning if scan is detected
    if length(ifits)>1
        fprintf('WARNING - guide_bot scan is not compatible with coatingwriter scans. This might result in errors!\n');
        pause(0.1);
    end
    
    % open ifit files to change price
    for j=1:length(ifits)
        fileID=fopen([char(ifits{j})],'r');
        tline = fgetl(fileID);
        k=1;
        A{k} = tline;
        while ischar(tline)
            k = k+1;
            tline = fgetl(fileID);
            A{k} = tline;
        end
        fclose(fileID);  
        
        % Insert new line in text array
        for k=1:length(A)
            if strfind(A{k},'scanname=')>0 
                if k<15
                    eval(A{k});
                    scanname=[scanname '_' num2str(thisPrice)]
                    A{k}=['scanname=''' scanname ''''];
                end
            end
            if strfind(A{k},'p.totalPrice=')>0
                A{k}=['p.totalPrice=''' num2str(thisPrice) '''; % total price (substrate plus coating)'];
            end
        end
        % Update file
        fileID=fopen([char(ifits{j})],'w');
        for j=1:length(A)
            fprintf(fileID,'%s\n',char(A{j}));
        end
        fclose(fileID);
    end
    


    %% Add line to cluster launch all file
    launch_cluster{end+1}=['cd '  deepestFolder '_price_' num2str(thisPrice)];
    launch_cluster{end+1}=['sbatch ' instrumentName '.batch'];
    launch_cluster{end+1}='cd ..';
    %launch_cluster{end+1}='clear';
    compile{end+1}=['cd ' deepestFolder '_price_' num2str(thisPrice)];
    compile{end+1}=['./compile_' instrumentName '.sh'];
    compile{end+1}='cd ..';
    compilePY{end+1}=['cd ' deepestFolder '_price_' num2str(thisPrice)];
    compilePY{end+1}=['./compile_' instrumentName '_py.sh'];
    compilePY{end+1}='cd ..';
    
    
    RunAll{end+1}=['cd ' deepestFolder '_price_' num2str(thisPrice)];
    RunAll{end+1}=['./compile_' instrumentName '.sh'];
    RunAll{end+1}=['sbatch ' instrumentName '.batch'];
    RunAll{end+1}='cd(thisFolder)';
    
    %% Add line to home-computer launch all file
    launch_home{end+1}=['% run number ' num2str(i) ];
    launch_home{end+1}=['cd ' deepestFolder '_price_' num2str(thisPrice)];
    for j=1:numel(ifits)
        launch_home{end+1}=['run(''' char(ifits{j}) ''')'];
    end
    launch_home{end+1}='cd(thisFolder)';
    
    %% Add line to home-computer launch all - parallel file
    if ParWorks==1
        launch_home_parallel{end+1}=['fileList{' num2str(i+1) '}=''' char(ifits{j}) ''';'];
        launch_home_parallel{end+1}=['pathList{' num2str(i+1) '}=[Folder slash ''' deepestFolder '_price_' num2str(thisPrice) '''];'];
    end
end

%% End launch files
launch_home_parallel{end+1}='';
launch_home_parallel{end+1}='parfor i=1:length(fileList)';
launch_home_parallel{end+1}='    parFun(pathList{i},fileList{i})';
launch_home_parallel{end+1}='end';
launch_home_parallel{end+1}='';
launch_home_parallel{end+1}=['fprintf(''All Done!\n'')'];

%% Print launch files
cd ..
% Cluster
fileID=fopen(['launch_' instrumentName '_scan.sh'],'w');
for i=1:length(launch_cluster)
    fprintf(fileID,'%s\n',char(launch_cluster{i}));
end
fclose(fileID);

fileID=fopen(['compile_all.sh'],'a+');
for i=1:length(compile)
    fprintf(fileID,'%s\n',char(compile{i}));
end
% Launch all:
fileID=fopen(['launch_all.sh'],'r+');
i = 1;
tline = fgetl(fileID);
launch_all_lines{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fileID);
    launch_all_lines{i} = tline;
end
launch_all_lines(end)=[];
fclose(fileID);
fileID=fopen(['launch_all.sh'],'w');
for i=1:length(launch_all_lines)
   if strcmp(launch_all_lines{i},['sbatch ./' instrumentName '.batch']) 
       launch_all_lines{i+1}=[''];
       launch_all_lines{i}=[''];
       launch_all_lines{i-1}=[''];
   end
   if strcmp(launch_all_lines{i},['watch -n 1 cat RunStatus.txt']) 
       launch_all_lines{i}=[''];
   end
end
for i=1:length(launch_all_lines)
    fprintf(fileID,'%s\n',char(launch_all_lines{i}));
end
fclose(fileID);

fileID=fopen(['launch_all.sh'],'a+');
for i=9:length(launch_cluster)
    fprintf(fileID,'%s\n',char(launch_cluster{i}));
end




fclose(fileID);
fileID=fopen(['compile_all_py.sh'],'a+');
for i=1:length(compilePY)
    fprintf(fileID,'%s\n',char(compilePY{i}));
end
fclose(fileID);

% Home
fileID=fopen(['launch_' instrumentName '_scan_home.m'],'w');
for i=1:length(launch_home)
    fprintf(fileID,'%s\n',char(launch_home{i}));
end
fclose(fileID);

% Home parallel
if ParWorks==1
    fileID=fopen(['launch_' instrumentName '_scan_home_parallel.m'],'w');
    for i=1:length(launch_home_parallel)
        fprintf(fileID,'%s\n',char(launch_home_parallel{i}));
    end
    fclose(fileID);
end
if ParWorks==1
if isfile(['launch_all_home_parallel.m'])
    fileID=fopen(['launch_all_home_parallel.m'],'a');
    for i=25:length(launch_home_parallel)
        fprintf(fileID,'%s\n',char(launch_home_parallel{i}));
    end
    fclose(fileID);
else
    fileID=fopen(['launch_all_home_parallel.m'],'w');
    for i=1:length(launch_home_parallel)
        fprintf(fileID,'%s\n',char(launch_home_parallel{i}));
    end
    fclose(fileID);
end
end

fprintf('All folders created, run using launch_%s.sh on cluster or launch_home.m on home computer\n',instrumentName)
end
