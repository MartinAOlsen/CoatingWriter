function InsertCoating(Coating,CoatingFirst,Declare,Define,Initialize,Monitor,lines_ifit,coatingOptions,Summary)
tic
runonce=0;
logOnce =0;
compileNopirce=0;
compileStepwise=0;
%CWfolder=cd;
CWfolder=mfilename('fullpath'); % Installation folder of CoatingWriter!
[upperPath, deepestFolder, ~] = fileparts(CWfolder);
CWfolder=upperPath;
[upperPath, deepestFolder, ~] = fileparts(coatingOptions.filePath);
thisFolder=upperPath;
% cd(CWfolder);
% cd PGE



%% Make path windows-compatible
nslash=length(strfind(coatingOptions.filePath,'/'));
nbackslash=length(strfind(coatingOptions.filePath,'\'));
if nslash<nbackslash
    slash='\';
else
    slash='/';
end

%% Start
fprintf('\nInserting CoatingWriter code in guide_bot files using InsertCoating.m:\n\n')

%% If +S is active, create some new instument files:
if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
    
    
end



%% Remove analysis if optimizer has build-in analysis:
if (strfind(coatingOptions.OptimizeMode,'+A')+strfind(coatingOptions.OptimizeMode,'+a'))>0
    coatingOptions.AnalyzeMode='none';
end


%% Shutter
if strcmp(coatingOptions.punishment,'shutter')
    shutter{1}='';
    shutter{2}='';
    shutter{3}='COMPONENT BS2 = Beamstop(xwidth=punishment/100000-1/100000, yheight=punishment/100000-1/100000, radius=0)';
    shutter{4}='AT (0,0, 0) RELATIVE PREVIOUS';
    shutter{5}='ROTATED (0,0,0) RELATIVE PREVIOUS';
    shutter{6}='';
    shutter{7}='COMPONENT EndOfelement_22= Arm()';
    shutter{8}='AT (0,0,u) RELATIVE PREVIOUS';
else
    shutter{1}='';
end




SepInd = findstr(coatingOptions.filePath, filesep);
nSep = length(SepInd);
InstrumentName = coatingOptions.filePath(SepInd(nSep)+1:length(coatingOptions.filePath));
folderPath = coatingOptions.filePath(1:SepInd(nSep));

P= mfilename('fullpath');
[P,~,~]=fileparts(P); 

copyfile([P slash 'guifiles' slash 'CoatingWriter_monitor.comp'],[coatingOptions.filePath slash 'CoatingWriter_monitor.comp']);
copyfile([P slash 'guifiles' slash 'CWGUI.m'],[coatingOptions.filePath slash 'CWGUI.m']);
copyfile([P slash 'guifiles' slash 'CW_optimizer.m'],[coatingOptions.filePath slash 'CW_optimizer.m']);
copyfile([P slash 'guifiles' slash 'LoadCWMonitor.m'],[coatingOptions.filePath slash 'LoadCWMonitor.m']);
copyfile([P slash 'guifiles' slash 'ManuallyChoosingOfOutput.m'],[coatingOptions.filePath slash 'ManuallyChoosingOfOutput.m']);
copyfile([P slash 'guifiles' slash 'optim_Genetic.m'],[coatingOptions.filePath slash 'optim_Genetic.m']);
copyfile([P slash 'guifiles' slash 'optim_Genetic_memory.m'],[coatingOptions.filePath slash 'optim_Genetic_memory.m']);
copyfile([P slash 'guifiles' slash 'optim_Genetic_memory_noWeights.m'],[coatingOptions.filePath slash 'optim_Genetic_memory_noWeights.m']);
copyfile([P slash 'guifiles' slash 'optim_Genetic_memory_prop.m'],[coatingOptions.filePath slash 'optim_Genetic_memory_prop.m']);
copyfile([P slash 'guifiles' slash 'optim_Genetic_weights.m'],[coatingOptions.filePath slash 'optim_Genetic_weights.m']);
copyfile([P slash 'guifiles' slash 'optim_random.m'],[coatingOptions.filePath slash 'optim_random.m']);
copyfile([P slash 'guifiles' slash 'optim_Walker.m'],[coatingOptions.filePath slash 'optim_Walker.m']);
copyfile([P slash 'guifiles' slash 'runMcStas.m'],[coatingOptions.filePath slash 'runMcStas.m']);


copyfile([P slash 'Insert_code' slash 'parFun.m'],[folderPath 'parFun.m']);
copyfile([P slash 'Insert_code' slash 'maxLine.m'],[folderPath 'maxLine.m']);
copyfile([P slash 'Insert_code' slash 'plot_guide.m'],[folderPath 'plot_guide.m']);
copyfile([P slash 'Insert_code' slash 'read_coating.m'],[folderPath 'read_coating.m']);
copyfile([P slash 'Insert_code' slash 'read_guide.m'],[folderPath 'read_guide.m']);


if strcmp(coatingOptions.scanType,'value') || strcmp(coatingOptions.scanType,'budget') || strcmp(coatingOptions.scanType,'pricescan') 
   copyfile([P slash 'Divergence_monitor_coating.comp'],[coatingOptions.filePath slash 'Divergence_monitor_coating.comp']);
   fprintf('Added Divergence_monitor_coating.comp for optimization of flux/price \n')
   coatingMonitor=1;
else
   coatingMonitor=0;
   try  % If the file is allready there, do nothing
       copyfile(['''' cd slash 'Divergence_monitor_coating.comp'''],['''' coatingOptions.filePath slash 'Divergence_monitor_coating.comp''']);
   end
end
parts = strsplit(coatingOptions.filePath, '/');
DirPart = parts{end};
if length(coatingOptions.filePath(1:end-length(DirPart)-1))<1
    copyfile([CWfolder slash 'Insert_code' slash 'AnalyseScan.m'],['AnalyseScan.m']);
    copyfile([CWfolder slash 'Insert_code' slash 'printStatusCW.m'],['printStatusCW.m']);
    copyfile([CWfolder slash 'Insert_code' slash 'RunStatus.txt'],['RunStatus.txt']);
else
    copyfile([CWfolder slash 'Insert_code' slash 'AnalyseScan.m'],[coatingOptions.filePath(1:end-length(DirPart)-1) slash 'AnalyseScan.m']);
    copyfile([CWfolder slash 'Insert_code' slash 'printStatusCW.m'],[coatingOptions.filePath slash 'printStatusCW.m']);
    copyfile([CWfolder slash 'Insert_code' slash 'RunStatus.txt'],[ coatingOptions.filePath(1:end-length(DirPart)-1) slash 'RunStatus.txt']);
end
%Insert Ifit

%% Find able instrumentfiles files in directory.
searchFileList{1}={'_optimize_ess.instr'};
searchFileList{2}={'_optimize.instr'};
searchFileList{3}={'_analyze_ess.instr'};
searchFileList{4}={'_analyze.instr'};

clear tmp;
for i=1:length(searchFileList)
    if exist([coatingOptions.filePath slash InstrumentName char(searchFileList{i})])>0
        tmp{i}=searchFileList{i};
    else
        tmp{i}=[];
    end
end
fileList = tmp(~cellfun(@isempty, tmp));

%% Add variables to the ifit file.
eval(sprintf('cd ''%s''',coatingOptions.filePath));
s = dir('*_ifit.m');
file_list = {s.name};
ifits={};
for i=1:length(file_list)
    tmp=char(file_list(i));
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
for k=1:length(ifits)
    if length(ifits{k})>0
        fileID=fopen([char(ifits{k})],'r');
        i = 1;
        tline = fgetl(fileID);
        A{i} = tline;
        while ischar(tline)
            i = i+1;
            tline = fgetl(fileID);
            A{i} = tline;
        end
        fclose(fileID);  
        
        fileID=fopen([CWfolder slash 'Insert_code' slash 'loadNSaveCW.m'],'r');  % add result saving lines to _ifit
        i = 1;
        tline = fgetl(fileID);
        writeResult{i} = tline;
        while ischar(tline)
            i = i+1;
            tline = fgetl(fileID);
            writeResult{i} = tline;
        end
        writeResult{end}='';
        fclose(fileID);  
        
        
        i=1;
        while i<length(A)-1
            tmp=A{i};
            if length(tmp)>2
                if strcmp('p.',tmp(1:2))
                    pLineLast=i;
                end
            end
            i=i+1;
        end
        before=A(1:pLineLast);
        after=A(pLineLast+1:length(A));
        
        % Insert all parts in a new list:
        B=[before,lines_ifit,after];
        
        
        %% FIX error in monitor output from new iFIT
        for i = 1:length(B) 
             if strfind(B{i},'optimal=monitor(1).Data.Parameters;')
                cutIndex = i;
            end 
        end
        clear tmp;
        tmp {1} = 'f = fieldnames(pars);'
        tmp {1} = 'optimal = p'
        tmp {3} = 'for i = 1:length(f)'
        tmp {4} = 'optimal.(f{i}) = pars.(f{i})'
        tmp {5} = 'end'
        
        B=[B(1:cutIndex-1),tmp,B(cutIndex+1:end)];
        
        % Split the "B" up into an "optimize" and "analyze"
        % part, and toggle the content from mode.
        
        % Find line numbers of interesting places
        i=1;
        while i<length(B)-1
            tmp=B{i};
            % The call to ifit:
            if strfind(tmp,['[pars,monitor,m,o]=mcstas'])
                pLineOptimize=i;
            end
            if strfind(tmp,['[pars_ideal,monitor_ideal,m_ideal,o_ideal]=mcstas'])
                pLineOptimize=i;
            end
            % The beginning of the analysis:
            if strfind(tmp,['%------------------------ Analy'])
                pLineAnalyzeBegin=i;
            end
            % The end of the analysis:
            if strfind(tmp,['fprintf(fid,[num2str(flux) '' = '' filename '' - '' inputstring ''\n''])'])
                pLineAnalyzeEnd=i+1;
            end
            i=i+1;
        end
        first=B(1:pLineOptimize-1);
        optimize=B(pLineOptimize);
        mid=B(pLineOptimize+1:pLineAnalyzeBegin-1);
        analyze=B(pLineAnalyzeBegin:pLineAnalyzeEnd)
        
        
        % Toggle Analyze modes
        if strcmp(coatingOptions.scanType,'manual')
            coatingOptions.AnalyzeMode = 'manual';
        end
        
        if strcmp(coatingOptions.AnalyzeMode,'standard') % Standard mode will run guide_bot analysis 
            
        elseif strcmp(coatingOptions.AnalyzeMode,'manual')
            button_analyze_script = analyze;
            analyze={};
            
            
        elseif strcmp(coatingOptions.AnalyzeMode,'reduced') || strcmp(coatingOptions.AnalyzeMode,'scan')   % Reduced mode will only save monitor_ALLW
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Analyze_Modes' slash 'reduced.m'],'r');  % add result saving lines to _ifit
            i = 1;
            analyze={}
            tline = fgetl(fileID);
            analyze{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                analyze{i} = tline;
            end
            fclose(fileID);  
            
        elseif strcmp(coatingOptions.AnalyzeMode,'none')  % Will remove the analysis    
            analyze={};
        end
        
        %% Toggle speedup mode
        if strcmp(coatingOptions.SpeedScanMode,'on')
            % Change neutron count:
            for i=1:length(first)
               if strfind(first{i},'options_cluster.ncount')
                    first{i}=['options_cluster.ncount=9*1e6;'];
               end
               if strfind(first{i},'options_single_cluster.ncount=')
                    first{i}=['options_single_cluster.ncount=1*1e7;'];
               end    
            end
            % Remove "exclusive" from batch file
            b = dir('*.batch');
            batch_list = {b.name};
            for i=1:length(batch_list)
                fileID=fopen([batch_list{i}],'r+');
                tline = fgetl(fileID);
                j=1;
                batch{j} = tline;   
                while ischar(tline)
                    j = j+1;
                    tline = fgetl(fileID);
                    if strfind(tline,'#SBATCH --exclusive');
                        batch{j} = '#SBATCH --exclusive';
                    else
                        batch{j} = tline;
                    end
%                     if strfind(tline,'NUMCORES=');
%                         batch{j} =['NUMCORES=2']
%                     end
                    if strfind(tline,'#SBATCH --job-name=');
                        batch{j} =['#SBATCH --job-name=' InstrumentName '_scan']
                    end
                    
                end
                fclose(fileID)
                fileID=fopen([batch_list{i}],'w');
                for j=1:length(batch)-1
                    if length(char(batch{j}))>0
                        fprintf(fileID,'%s\n',char(batch{j}));
                    end
                end
                fclose(fileID);
            end
        end
            
        %% Make a list of coating-parameters
        coatingVars=sprintf('coatingVars={');
        for i=1:length(lines_ifit)
            if strfind(lines_ifit{i},'p.')
                tempLines=strsplit(lines_ifit{i},'=');
                tempLines=tempLines{1};
                coatingVars=sprintf('%s''%s'',',coatingVars,tempLines(3:end));
                
            end

          
        end
        coatingVars=sprintf('%s};',coatingVars(1:end-1));


	if strcmp(coatingOptions.scanType,'manual')
		coatingOptions.OptimizeMode = 'manual';
	end
        %% Toggle Optimization modes
        if strfind(coatingOptions.OptimizeMode,'standard') % Standard mode will run guide_bot analysis
            % Load script:
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'standard.m'],'r');  % add result saving lines to _ifit
            i = 1;
            tline = fgetl(fileID);
            optimize{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                optimize{i} = tline;
            end
            fclose(fileID); 
            % Insert coating variable list:
            optimize=[['mode=''' coatingOptions.scanType ''';'],optimize];
            i=1;
            
            while i < length(optimize)
                tmp=optimize{i};
                % Insert stepwise scan code:
                if strfind(tmp,['%% ADDSTEPSCAN %%']) 
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'stepwiseScan.m'],'r');  % add result saving lines to _ifit
                            linenr = 1;
                            tline = fgetl(fileID);
                            addStepScan{linenr} = tline;
                            while ischar(tline)
                                linenr = linenr+1;
                                tline = fgetl(fileID);
                                addStepScan{linenr} = tline;
                            end
                            addStepScan{end}=[];
                            fclose(fileID);
                        optimize=[optimize(1:i-1) addStepScan optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];                   
                    end
                end
                
                i=i+1;
            end
            
        elseif strfind(coatingOptions.OptimizeMode,'manual') % Standard mode will run guide_bot analysis
            % Load script:
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'manual.m'],'r');  % add result saving lines to _ifit
            i = 1;
            tline = fgetl(fileID);
            optimize{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                optimize{i} = tline;
            end
            fclose(fileID); 
            % Insert coating variable list:
            optimize=[['mode=''' coatingOptions.scanType ''';'],optimize];
            i=1;
            
        
           
            
        elseif strfind(coatingOptions.OptimizeMode,'3step')
            % Load script:
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash '3step.m'],'r');  % add result saving lines to _ifit
            i = 1;
            tline = fgetl(fileID);
            optimize{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                optimize{i} = tline;
            end
            optimize{end}=[];
            fclose(fileID); 
            % Insert coating variable list:
            optimize=[['mode=''' coatingOptions.scanType ''';'],optimize];
            i=1;
            
            while i < length(optimize)
                tmp=optimize{i};
                %% Insert coating variables
                if strfind(tmp,['% INSERT COATING VARIABLES HERE!'])
                    optimize{i}=coatingVars;
                end
                %% Insert first instrument _noPrice if +C mode is active
                if strfind(tmp,['% INSERT FIRST INSTRUMENT HERE!'])
                    if sum([strfind(coatingOptions.OptimizeMode,'+C'),strfind(coatingOptions.OptimizeMode,'+c')])>0
                        optimize{i-1}=['mode1=''coating'';']
                        optimize{i}=['[pars,monitor,m,o]=mcstas([instrument_name ''_optimize_noPrice.instr''],tempPar,options{select});']
                    else
                        optimize{i-1}=['mode1=''value'';']
                        optimize{i}=['[pars,monitor,m,o]=mcstas([instrument_name ''_optimize.instr''],tempPar,options{select});']
                    end
                end
                %% Insert analysis if +A mode is active
                if strfind(tmp,['%% ADDANALYSIS %%'])
                    if sum([strfind(coatingOptions.OptimizeMode,'+A'),strfind(coatingOptions.OptimizeMode,'+a')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'addAnalysis.m'],'r');  % add result saving lines to _ifit
                        linenr = 1;
                        tline = fgetl(fileID);
                        addAnalysis{linenr} = tline;
                        while ischar(tline)
                            linenr = linenr+1;
                            tline = fgetl(fileID);
                            addAnalysis{linenr} = tline;
                        end
                        addAnalysis{end}=[];
                        fclose(fileID);
                        optimize=[optimize(1:i-1) addAnalysis optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];
                    end
                end
                %% Insert stepScan if +S mode is active
                % Add one to list of steps
                if strfind(tmp,['steps='])>0 %&& strfind(tmp,['+'])==0
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        optimize{i}=[tmp(1:end-1) '+1;'];
                    end
                end
                % Insert stepwise scan code:
                if strfind(tmp,['%% ADDSTEPSCAN %%']) 
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'stepwiseScan.m'],'r');  % add result saving lines to _ifit
                            linenr = 1;
                            tline = fgetl(fileID);
                            addStepScan{linenr} = tline;
                            while ischar(tline)
                                linenr = linenr+1;
                                tline = fgetl(fileID);
                                addStepScan{linenr} = tline;
                            end
                            addStepScan{end}=[];
                            fclose(fileID);
                        optimize=[optimize(1:i-1) addStepScan optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];                   
                    end
                end
               
                i=i+1;
            end
            
            
        elseif strfind(coatingOptions.OptimizeMode,'6step')
            % Load script:
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash '6step.m'],'r');  % add result saving lines to _ifit
            i = 1;
            tline = fgetl(fileID);
            optimize{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                optimize{i} = tline;
            end
            optimize{end}=[];
            fclose(fileID); 
            % Insert coating variable list:
            i=1;
            optimize=[['mode=''' coatingOptions.scanType ''';'],optimize];
            while i < length(optimize)
                tmp=optimize{i};
                %% Insert coating variables
                if strfind(tmp,['% INSERT COATING VARIABLES HERE!'])
                    optimize{i}=coatingVars;
                end
               %% Insert first instrument _noPrice if +C mode is active
                if strfind(tmp,['% INSERT FIRST INSTRUMENT HERE!'])
                    if sum([strfind(coatingOptions.OptimizeMode,'+C'),strfind(coatingOptions.OptimizeMode,'+c')])>0
                        optimize{i-1}=['mode1=''coating'';']
                        optimize{i}=['[pars,monitor,m,o]=mcstas([instrument_name ''_optimize_noPrice.instr''],tempPar,options{select});']
                    else
                        optimize{i-1}=['mode1=''value'';']
                        optimize{i}=['[pars,monitor,m,o]=mcstas([instrument_name ''_optimize.instr''],tempPar,options{select});']
                    end
                end
                %% Insert analysis if +A mode is active
                if strfind(tmp,['%% ADDANALYSIS %%'])
                    if sum([strfind(coatingOptions.OptimizeMode,'+A'),strfind(coatingOptions.OptimizeMode,'+a')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'addAnalysis.m'],'r');  % add result saving lines to _ifit
                        linenr = 1;
                        tline = fgetl(fileID);
                        addAnalysis{linenr} = tline;
                        while ischar(tline)
                            linenr = linenr+1;
                            tline = fgetl(fileID);
                            addAnalysis{linenr} = tline;
                        end
                        addAnalysis{end}=[];
                        fclose(fileID);
                        optimize=[optimize(1:i-1) addAnalysis optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];
                    end
                end



		%% Insert stepScan if +S mode is active
                % Add one to list of steps
                if strfind(tmp,['steps='])>0 %&& strfind(tmp,['+'])==0
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        optimize{i}=[tmp(1:end-1) '+1;'];
                    end
                end
                % Insert stepwise scan code:
                if strfind(tmp,['%% ADDSTEPSCAN %%']) 
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'stepwiseScan.m'],'r');  % add result saving lines to _ifit
                            linenr = 1;
                            tline = fgetl(fileID);
                            addStepScan{linenr} = tline;
                            while ischar(tline)
                                linenr = linenr+1;
                                tline = fgetl(fileID);
                                addStepScan{linenr} = tline;
                            end
                            addStepScan{end}=[];
                            fclose(fileID);
                        optimize=[optimize(1:i-1) addStepScan optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];                   
                    end
                end
                i=i+1;
            end    
            
        elseif strfind(coatingOptions.OptimizeMode,'9step')
            % Load script:
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash '9step.m'],'r');  % add result saving lines to _ifit
            i = 1;
            tline = fgetl(fileID);
            optimize{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                optimize{i} = tline;
            end
            optimize{end}=[];
            fclose(fileID); 
            % Insert coating variable list:
            i=1;
            optimize=[['mode=''' coatingOptions.scanType ''';'],optimize];
            while i < length(optimize)
                tmp=optimize{i};
                %% Insert coating variables
                if strfind(tmp,['% INSERT COATING VARIABLES HERE!'])
                    optimize{i}=coatingVars;
                end
                %% Insert first instrument _noPrice if +C mode is active
                if strfind(tmp,['% INSERT FIRST INSTRUMENT HERE!'])
                    if sum([strfind(coatingOptions.OptimizeMode,'+C'),strfind(coatingOptions.OptimizeMode,'+c')])>0
                        optimize{i-1}=['mode1=''coating'';']
                        optimize{i}=['[pars,monitor,m,o]=mcstas([instrument_name ''_optimize_noPrice.instr''],tempPar,options{select});']
                    else
                        optimize{i-1}=['mode1=''value'';']
                        optimize{i}=['[pars,monitor,m,o]=mcstas([instrument_name ''_optimize.instr''],tempPar,options{select});']
                    end
                end
                %% Insert analysis if +A mode is active
                if strfind(tmp,['%% ADDANALYSIS %%'])
                    if sum([strfind(coatingOptions.OptimizeMode,'+A'),strfind(coatingOptions.OptimizeMode,'+a')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'addAnalysis.m'],'r');  % add result saving lines to _ifit
                        linenr = 1;
                        tline = fgetl(fileID);
                        addAnalysis{linenr} = tline;
                        while ischar(tline)
                            linenr = linenr+1;
                            tline = fgetl(fileID);
                            addAnalysis{linenr} = tline;
                        end
                        addAnalysis{end}=[];
                        fclose(fileID);
                        optimize=[optimize(1:i-1) addAnalysis optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];
                    end
                end
%                 %% Add print-logic
%                 if strfind(tmp,'printStatusCW')
%                     optionsIndex=strfind(tmp,'printStatusCW')+length('printStatusCW');
%                     options=eval(tmp(optionsIndex:end))
%                     optimize=[optimize(1:i-1) addPrint optimize(i+1:end)];
%                 end
                %% Insert stepScan if +S mode is active
                % Add one to list of steps
                if strfind(tmp,['steps='])>0 %&& strfind(tmp,['+'])==0
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        optimize{i}=[tmp(1:end-1) '+1;'];
                    end
                end
                % Insert stepwise scan code:
                if strfind(tmp,['%% ADDSTEPSCAN %%']) 
                    if sum([strfind(coatingOptions.OptimizeMode,'+S'),strfind(coatingOptions.OptimizeMode,'+s')])>0
                        fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'stepwiseScan.m'],'r');  % add result saving lines to _ifit
                            linenr = 1;
                            tline = fgetl(fileID);
                            addStepScan{linenr} = tline;
                            while ischar(tline)
                                linenr = linenr+1;
                                tline = fgetl(fileID);
                                addStepScan{linenr} = tline;
                            end
                            addStepScan{end}=[];
                            fclose(fileID);
                        optimize=[optimize(1:i-1) addStepScan optimize(i+1:end)];
                    else
                        optimize=[optimize(1:i-1) optimize(i+1:end)];                   
                    end
                end
                
                i=i+1;
            end
            
            
        elseif strcmp(coatingOptions.OptimizeMode,'singleParScan')
            % Load script:
            fileID=fopen([CWfolder slash 'Insert_code' slash 'Optimize_Modes' slash 'singleParScan.m'],'r');  % add result saving lines to _ifit
            i = 1;
            tline = fgetl(fileID);
            optimize{i} = tline;
            while ischar(tline)
                i = i+1;
                tline = fgetl(fileID);
                optimize{i} = tline;
            end
            fclose(fileID); 
        end
        
        if strcmp(coatingOptions.scanType,'manual') == 0
            tryL{1}=['try'];
        else
            tryL{1}=[''];
        end
        
        printFinish{1}=['%%%%%%%%%% Print'];
        printFinish{2}=['printStep=''finish'';'];
        printFinish{3}=['time=toc;'];
        printFinish{4}=['eval(printScript)'];
        printFinish{5}=['%%%%%%%%%%'];
        if strcmp(coatingOptions.scanType,'manual')==0
        printFinish{6}=['catch ME'];
        printFinish{7}=['    %%%%%%%%%% Print'];
        printFinish{8}=['    printStep=''failed'';'];
        printFinish{9}=['    time=toc;'];
        printFinish{10}=['    eval(printScript)'];
        printFinish{11}=['    %%%%%%%%%%'];
        printFinish{12}=['    rethrow(ME)'];
        printFinish{13}=['end'];
        end
         
        
        if sum([sum(coatingOptions.AnalyzeStepwise==1),strcmp(coatingOptions.AnalyzeStepwise,'on')]) >0    %% Stepwise 
            stepwise{1}=['options_single{2}.ncount=5*1e7'];

            stepwise{end+1}=['    optimal.WaveMin=p.WaveMin;'];
            stepwise{end+1}=['    optimal.WaveMax=p.WaveMax;'];
            stepwise{end+1}=['    options_single{select}.dir=[filename ''stepwise''];']
            stepwise{end+1}='    optimal_visualizer = optimal;';
            stepwise{end+1}=['    optimal_visualizer.file_name = [filename ''_geometry.dat''];'];
            stepwise{end+1}=['    monitor_stepwise=mcstas([instrument_name ''_analyze_stepwise.instr''],optimal_visualizer,options_single{select});'];
            stepwise{end+1}=['    Result_this.stepMon=monitor_stepwise;']; 

            stepwise{end+1}=['    optimal.WaveMin=0.45;'];
            stepwise{end+1}=['    optimal.WaveMax=0.45;'];
            stepwise{end+1}=['    options_single{select}.dir=[filename ''stepwise_epithermal''];']
            stepwise{end+1}='    optimal_visualizer = optimal;';
            stepwise{end+1}=['    optimal_visualizer.file_name = [filename ''_geometry.dat''];'];
            stepwise{end+1}=['    monitor_stepwise=mcstas([instrument_name ''_analyze_stepwise.instr''],optimal_visualizer,options_single{select});'];
            stepwise{end+1}=['    Result_this.stepMon_epithermal=monitor_stepwise;'];
            
            stepwise{end+1}=['    optimal.WaveMin=2;'];
            stepwise{end+1}=['    optimal.WaveMax=2;'];
            stepwise{end+1}=['    options_single{select}.dir=[filename ''stepwise_thermal''];']
            stepwise{end+1}='    optimal_visualizer = optimal;';
            stepwise{end+1}=['    optimal_visualizer.file_name = [filename ''_geometry.dat''];'];
            stepwise{end+1}=['    monitor_stepwise=mcstas([instrument_name ''_analyze_stepwise.instr''],optimal_visualizer,options_single{select});'];
            stepwise{end+1}=['    Result_this.stepMon_thermal=monitor_stepwise;'];
            
            stepwise{end+1}=['    optimal.WaveMin=4;'];
            stepwise{end+1}=['    optimal.WaveMax=4;'];
            stepwise{end+1}=['    options_single{select}.dir=[filename ''stepwise_cold''];']
            stepwise{end+1}='    optimal_visualizer = optimal;';
            stepwise{end+1}=['    optimal_visualizer.file_name = [filename ''_geometry.dat''];'];
            stepwise{end+1}=['    monitor_stepwise=mcstas([instrument_name ''_analyze_stepwise.instr''],optimal_visualizer,options_single{select});'];
            stepwise{end+1}=['    Result_this.stepMon_cold=monitor_stepwise;'];
            

            newIfitFile=[tryL,first,optimize,mid,analyze,stepwise,printFinish,writeResult]
        else
            newIfitFile=[tryL,first,optimize,mid,analyze,printFinish,writeResult]
        end
        

        %fileID=fopen([coatingOptions.filePath '/' char(ifits(k))],'w');
        fileID=fopen([char(ifits{k})],'w');
        for i=1:length(newIfitFile)
            if isa(newIfitFile{i},'char')==0
                newIfitFile(i)={''}
            end
            fprintf(fileID,'%s\n',char(newIfitFile{i}));
        end
        fclose(fileID);
        fprintf('Inserted %i lines in: %s\n',length(lines_ifit),char(ifits{k}))
    end
end

%% The lines that read the correction array if stepScan is activated
[loadStepScan] = StepScanLoader(coatingOptions);



%% Fix the instrument files
for instrFile=1:length(fileList)
    isOptimization=0;
    if sum(strfind(char(fileList{instrFile}),'Optimize'))>0 || sum(strfind(char(fileList{instrFile}),'optimize'))>0
         isOptimization=1;
    end
    clear A
    fileID=fopen([coatingOptions.filePath slash InstrumentName char(fileList{instrFile})],'r');
    i = 1;
    tline = fgetl(fileID);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fileID);
        A{i} = tline;
    end
    fclose(fileID);
    i=1;
    for j=1:coatingOptions.nInput
        localOptions=eval(['coatingOptions.Input' num2str(j)]);
        tmpSearch=['m=m' num2str(j) ','];

        tmpSearchCurveA=['ma=m' num2str(j)];
        tmpSearchCurveI=['mi=m' num2str(j)];
        tmpSearchCurveS=['ms=m' num2str(j)];
        i=1;
        %% Loop through all lines in .instr file, and insert corrections.
        while i<length(A)-1
            if strcmp(coatingOptions.scanType,'manual') &&  sum(strfind(char(fileList{instrFile}),'optimize'))>0
                if strfind(A{i},'COMPONENT Div2d_sample_B = Divergence_monitor(')>0
                    A{i}= ['COMPONENT CW_mon = CoatingWriter_monitor(filename = "Monitor.txt", xwidth = sizeX,yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y,Price=sectionPrices+TotalSubstratePrice,min_lambda=1.5,max_lambda=WaveMax,background_min_lambda=' num2str(coatingOptions.WavelengthBackgroundMin)  ',background_max_lambda=' num2str(coatingOptions.WavelengthBackgroundMax) ')'];
                end
                if strfind(A{i},'    nh = 20, nv = 20, filename = "Div2d_sample_B", xwidth = sizeX,')>0
                    A{i}= 'AT (0, 0, sample_dist) RELATIVE PREVIOUS';
                end
                if strfind(A{i},'yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y)')>0
                     A{i}='';
                end
                if strfind(A{i},'AT (0, 0,sample_dist) RELATIVE PREVIOUS')>0
                     A{i}=''
                end
                if strfind(A{i},'dLambda = 0.5*(WaveMax - WaveMin);')> 0
                    A{i} = ['dLambda = 0.5*(WaveMax - ' num2str(coatingOptions.WavelengthBackgroundMin) ');']
                end
                if strfind(A{i},'Lambda0 = dLambda+WaveMin;')> 0
                    if coatingOptions.WavelengthBackgroundMin == 0
                        A{i} = ['Lambda0 = dLambda+' num2str(coatingOptions.WavelengthBackgroundMin+0.00001) ';']
                    else
                        A{i} = ['Lambda0 = dLambda+' num2str(coatingOptions.WavelengthBackgroundMin) ';']
                    end
                end
            end
            
            
            
           if strfind(A{i},tmpSearch)>0
               tmp=A{i};
               before=tmp(1:strfind(A{i},tmpSearch)-1);
               removeLength=length(num2str(j))+3;
               after=tmp(strfind(A{i},tmpSearch)+removeLength:length(tmp));


               lock=eval(sprintf('coatingOptions.Input%i.fixSides',j));
               % Correct names of the coating variables of the most common
               % guide geometries (P,S and E input). Dependant on what
               % sides are locked together, look at variable "lock"
               if strcmp(lock,'HV')
                   if strcmp(localOptions.type,'S')
                        insert=[before 'mright = mValues' num2str(j) 'horizontal, mleft = mValues' num2str(j) 'horizontal,mtop = mValues' num2str(j) 'vertical,mbottom = mValues' num2str(j) 'vertical' after];
                   else
                        insert=[before 'mvaluesright = mValues' num2str(j) 'horizontal, mvaluesleft = mValues' num2str(j) 'horizontal,mvaluestop = mValues' num2str(j) 'vertical,mvaluesbottom = mValues' num2str(j) 'vertical,seglength = elementLength' num2str(j) after];
                   end
               end
               if strcmp(lock,'none')
                   if strcmp(localOptions.type,'S')
                        insert=[before 'mright = mValues' num2str(j) 'right, mleft = mValues' num2str(j) 'left,mtop = mValues' num2str(j) 'top,mbottom = mValues' num2str(j) 'bottom' after];
                   else
                        insert=[before 'mvaluesright = mValues' num2str(j) 'right, mvaluesleft = mValues' num2str(j) 'left,mvaluestop = mValues' num2str(j) 'top,mvaluesbottom = mValues' num2str(j) 'bottom,seglength = elementLength' num2str(j) after];
                   end
               end
               if strcmp(lock,'all')
                   if strcmp(localOptions.type,'S')
                        insert=[before 'mright = mValues' num2str(j) ', mleft = mValues' num2str(j) ',mtop = mValues' num2str(j) ',mbottom = mValues' num2str(j) after];
                   else
                        insert=[before 'mvaluesright = mValues' num2str(j) ', mvaluesleft = mValues' num2str(j) ',mvaluestop = mValues' num2str(j) ',mvaluesbottom = mValues' num2str(j) ',seglength = elementLength' num2str(j) after];
                   end
               end
               if strcmp(lock,'curve')
                   if strcmp(localOptions.type,'S')
                        insert=[before 'mright = mValues' num2str(j) 'right, mleft = mValues' num2str(j) 'left,mtop = mValues' num2str(j) 'vertical,mbottom = mValues' num2str(j) 'vertical' after];
                   else
                        insert=[before 'mvaluesright = mValues' num2str(j) 'right, mvaluesleft = mValues' num2str(j) 'left,mvaluestop = mValues' num2str(j) 'vertical,mvaluesbottom = mValues' num2str(j) 'vertical,seglength = elementLength' num2str(j) after];
                   end
               end
               A{i}=insert;
           end
           % Insert the correct names for the mvalues of the curved
           % sections (mi, ma, ms): 
           if strfind(A{i},tmpSearchCurveA)>0
               tmp=A{i};
               before=tmp(1:strfind(A{i},tmpSearchCurveA)-1);
               after=tmp(strfind(A{i},tmpSearchCurveA)+5:length(tmp));
               lock=eval(sprintf('coatingOptions.Input%i.fixSides',j));
               insert=[before 'ma=mValues' num2str(j) 'outside' after];
               A{i}=insert;
           end
            if strfind(A{i},tmpSearchCurveI)>0
               tmp=A{i};
               before=tmp(1:strfind(A{i},tmpSearchCurveI)-1);
               after=tmp(strfind(A{i},tmpSearchCurveI)+5:length(tmp));
               lock=eval(sprintf('coatingOptions.Input%i.fixSides',j));
               insert=[before 'mi=mValues' num2str(j) 'inside' after];
               A{i}=insert;
            end
            if strfind(A{i},tmpSearchCurveS)>0
               tmp=A{i};
               before=tmp(1:strfind(A{i},tmpSearchCurveS)-1);
               after=tmp(strfind(A{i},tmpSearchCurveS)+5:length(tmp));
               lock=eval(sprintf('coatingOptions.Input%i.fixSides',j));
               insert=[before 'ms=mValues' num2str(j) 'vertical' after];
               A{i}=insert;
            end

            % Make the optimizer call the monitor that includes price
            % instead of the regular divergence monitor, if needed:
            if isOptimization==1 && coatingMonitor==1; 
                if sum(strfind(A{i},'Divergence_monitor('))>0        
                   tmp=A{i};
                   before=tmp(1:strfind(A{i},'Divergence_monitor')-1);
                   after=tmp(strfind(A{i},'Divergence_monitor')+18:length(tmp));
                   insert=[before 'Divergence_monitor_coating' after];
                   A{i}=insert;
                end

                 if sum(strfind(A{i},'maxdiv_v = divreq_y)'))>0        
                   tmp=A{i};
                   before=tmp(1:strfind(A{i},'maxdiv_v = divreq_y)')+18);
                   after=tmp(strfind(A{i},'maxdiv_v = divreq_y)')+20:length(tmp));
                   insert=[before ',Price=sectionPrices+TotalSubstratePrice)' after];
                   A{i}=insert;
                end
            end
           i=i+1; 
        end
    end

    i=1;
    while i<length(A)-1
        tmp=A{i};
        if length(tmp)>16
            if strcmp('DEFINE INSTRUMENT',tmp(1:17));
                firstBreakAfter=i;
            end
        end
        if length(tmp)>6
            if strcmp('DECLARE',tmp(1:7));
                secondBreakAfter=i+1;
            end
        end
        if length(tmp)>9
            if strcmp('INITIALIZE',tmp(1:10));
                thirdBreakAfter=i+1;
            end
        end
        if length(tmp)>4
            if strcmp('TRACE',tmp(1:5));
                fourthBreakAfter=i-2;
            end
        end
        if length(tmp)>35
            if strcmp('AT (0,0,guide_start) RELATIVE source',tmp(1:36));
                    shutterBreakAfter=i;
            end
        end
        i=i+1;
    end
    if exist('shutterBreakAfter')==0
        shutterBreakAfter=fourthBreakAfter+2;
    end

    if strcmp(coatingOptions.scanType,'budget') || strcmp(coatingOptions.scanType,'pricescan')
        if isOptimization==1
            firstSnip=A(1:firstBreakAfter);
            secondSnip=A(firstBreakAfter+1:secondBreakAfter);
            thirdSnip=A(secondBreakAfter+1:thirdBreakAfter);
            fourthSnip=A(thirdBreakAfter+1:fourthBreakAfter);
            fifthSnip=A(fourthBreakAfter+1:shutterBreakAfter);
            sixthSnip=A(shutterBreakAfter+1:length(A));
           
            
        else
            firstSnip=A(1:firstBreakAfter);
            secondSnip=A(firstBreakAfter+1:secondBreakAfter);
            thirdSnip=A(secondBreakAfter+1:thirdBreakAfter);
            fourthSnip=A(thirdBreakAfter+1:fourthBreakAfter);
            sixthSnip=A(fourthBreakAfter+1:length(A));
        end
        
        %Quick Fix
        if strcmp(sixthSnip{end},'END')==0
           sixthSnip{end}={} ;
        end
        
        if isOptimization==1
            %% Make the optimizer save the price each run, for making an
            % envelope plot:
            if logOnce == 0
                out1{1}='char fileup[100];';
                out1{end+1}='char fileP[100];';
                out1{end+1}='sprintf(fileup,"%s%s%s","priceList",scanname,".txt");';
                out1{end+1}='sprintf(fileP,"%s%s%s","priceListPunished",scanname,".txt");';
                out1{end+1}='MPI_MASTER(';
                out1{end+1}=['FILE *priceFile=fopen(fileup,"a");'];
                out1{end+1}=['fprintf(priceFile,"%2.5f\n",sectionPricesUnpunished);'];
                out1{end+1}='fclose(priceFile);';
                out1{end+1}=['FILE *priceFilePunished=fopen(fileP,"a");'];
                out1{end+1}=['fprintf(priceFilePunished,"%2.5f\n",sectionPrices);'];
                out1{end+1}='fclose(priceFile);';
                out1{end+1}=')';
                logOnce=1;
                
% 
%                 out2{1}='MPI_MASTER(';
%                 out2{end+1}=['FILE *priceFile=fopen("priceList.txt","a");'];
%                 out2{end+1}=['fprintf(priceFile,"%2.5f\n",sectionPrices);'];
%                 out2{end+1}='fclose(priceFile);';
%                 out2{end+1}=['FILE *priceFilePunished=fopen("priceListPunished.txt","a");'];
%                 out2{end+1}=['fprintf(priceFilePunished,"%2.5f\n",sectionPrices);'];
%                 out2{end+1}='fclose(priceFile);';
%                 out2{end+1}=')';
            end
   
            %% Add punishment
            if strcmp(coatingOptions.punishment,'potential') && runonce==0;
               Punishment{1} =['double Error=fabs(sectionPrices-totalPrice);'];
               Punishment{end+1} =['double sectionPricesUnpunished =sectionPrices;']
               Punishment{end+1} =['if (Error>' num2str(coatingOptions.budgetMaxError) ') {']
               Punishment{end+1} =['sectionPrices+=pow((Error-' num2str(coatingOptions.budgetMaxError) '),' num2str(coatingOptions.potentialPower) ');'];
               Punishment{end+1} ='}'
               if strcmp(coatingOptions.SpeedScanMode,'on')
                   Punishment{end+1} =['double orgNCount = mcget_ncount();'];
                   Punishment{end+1} =['int newNCount = 1e4+round((orgNCount)*exp(-(Error/200)));'];
                   Punishment{end+1} =['mcset_ncount(newNCount);'];
                   Punishment{end+1} =['printf("-----NEUTRONS=%i (changed %2.2f %% from %2.2f)----\n",newNCount,-100*(((orgNCount-newNCount)/orgNCount)),orgNCount);'];
               end
               runonce=1;
            end
            
            newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,loadStepScan,fourthSnip,CoatingFirst,Coating,Punishment,out1,fifthSnip,shutter,sixthSnip];
        else
            if sum(strfind(char(fileList{instrFile}),'ess'))>0
                newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,loadStepScan,fourthSnip,CoatingFirst,Coating,sixthSnip];
            else
                newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,loadStepScan,fourthSnip,CoatingFirst,Coating,Summary,sixthSnip];
            end
        end
        
    else
        firstSnip=A(1:firstBreakAfter);
        secondSnip=A(firstBreakAfter+1:secondBreakAfter);
        thirdSnip=A(secondBreakAfter+1:thirdBreakAfter);
        fourthSnip=A(thirdBreakAfter+1:fourthBreakAfter);
        fifthSnip=A(fourthBreakAfter+1:length(A));
        
        
        if strcmp(fifthSnip{end},'END')==0
           fifthSnip{end}={} ;
        end
        
        
        if isOptimization==1
            %% Make the optimizer save the price each run, for making an
            % envelope plot:
            if logOnce == 0
                out1{1}='char fileup[100];';
                out1{end+1}='char fileP[100];';
                out1{end+1}='sprintf(fileup,"%s%s%s","priceList",scanname,".txt");';
                out1{end+1}='sprintf(fileP,"%s%s%s","priceListPunished",scanname,".txt");';
                out1{end+1}='MPI_MASTER(';
                out1{end+1}=['FILE *priceFile=fopen(fileup,"a");'];
                out1{end+1}=['fprintf(priceFile,"%2.5f\n",sectionPrices);'];
                out1{end+1}='fclose(priceFile);';
                out1{end+1}=['FILE *priceFilePunished=fopen(fileP,"a");'];
                out1{end+1}=['fprintf(priceFilePunished,"%2.5f\n",sectionPrices);'];
                out1{end+1}='fclose(priceFile);';
%                 out1{end+1}=')';
                logOnce=1;
            end
            
            newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,loadStepScan,fourthSnip,CoatingFirst,Coating,out1,Summary,')',fifthSnip];
        else
            if sum(strfind(char(fileList{instrFile}),'ess'))>0
                newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,loadStepScan,fourthSnip,CoatingFirst,Coating,fifthSnip];
            else
                 newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,loadStepScan,fourthSnip,CoatingFirst,Coating,Summary,fifthSnip];
               %newOptimFile=[firstSnip,Define,secondSnip,Declare,thirdSnip,Initialize,fourthSnip,CoatingFirst,Coating,fifthSnip];
            end
        end

    end
%     if isOptimization==1
% 
%     end




    

    fileID=fopen([coatingOptions.filePath '/' InstrumentName char(fileList{instrFile})],'w');
    for i=1:length(newOptimFile)
        fprintf(fileID,'%s\n',char(newOptimFile{i}));
    end
    fclose(fileID);
    fprintf('Inserted %i lines in: %s\n',length(Define)+length(Declare)+length(Initialize)+length(CoatingFirst)+length(Coating),[InstrumentName char(fileList{instrFile})])

end


if sum(coatingOptions.AnalyzeStepwise) > 0 || strcmp(coatingOptions.AnalyzeStepwise,'on')   %% Stepwise 
        % Dublicate analyze file
        copyfile([folderPath '/' InstrumentName '/' InstrumentName '_analyze.instr'],[folderPath '/' InstrumentName '/' InstrumentName '_analyze_stepwise.instr'])
        
        % Load file
        fileID=fopen([folderPath '/' InstrumentName '/' InstrumentName '_analyze_stepwise.instr'],'r');
        i = 1;
        tline = fgetl(fileID);
        A{i} = tline;
        while ischar(tline)
            i = i+1;
            tline = fgetl(fileID);
            A{i} = tline;
        end
        A{end}=[];
        fclose(fileID);
        
        % Remove old monitors
        for i=1:length(A)
           if strfind(A{i},'COMPONENT Lmon_guide_end = L_monitor(') 
                StartLine=i;
           end
           if strfind(A{i},'FINALLY') 
                EndLine=i-1;
           end
        end
        A(StartLine:EndLine)=[];
        % Add new monitors
        monitorPlacementList=[];
        for i=1:length(A)
             if sum(strfind(A{i},'COMPONENT EndOfelement'))>0 || sum(strfind(A{i},'COMPONENT StartOfGuide'))>0 
                 monitorPlacementList(end+1)=i+3;
             end
        end
        
        for i=1:length(monitorPlacementList)
            clear monLine;
            % Intensity:
            monLine{1}=['COMPONENT Lmon_sample_' num2str(i) ' = L_monitor('];
            monLine{end+1}=['nL = 1000, filename = "Lmon_sample_' num2str(i) '", xwidth = 1.0, restore_neutron=1,yheight = 1.0, Lmin = 0, Lmax = 10)'];
            monLine{end+1}=['AT (0, 0,u) RELATIVE PREVIOUS'];
            monLine{end+1}='';
            
            % Brilliance transfer
            monLine{end+1}=['COMPONENT Div2d_sample_B_' num2str(i) ' = Divergence_monitor('];
            if i<length(monitorPlacementList)
                monLine{end+1}=['nh = 200, nv = 200, filename = "Div2d_sample_B_' num2str(i) '", xwidth = startx' num2str(length(InstrumentName)-(i-1)) ','];
                monLine{end+1}=['yheight = startx' num2str(length(InstrumentName)-(i-1)) ', maxdiv_h = divreq_x*(sizeX/startx' num2str(length(InstrumentName)-(i-1)) '), maxdiv_v = divreq_y*(sizeY/starty' num2str(length(InstrumentName)-(i-1)) '),restore_neutron=1)'];
            else
                monLine{end+1}=['nh = 200, nv = 200, filename = "Div2d_sample_B_' num2str(i) '", xwidth = sizeX,'];
                monLine{end+1}=['yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y,restore_neutron=1)'];
            end
            monLine{end+1}=['AT (0, 0,sample_dist) RELATIVE PREVIOUS'];
            monLine{end+1}=['EXTEND'];
            monLine{end+1}=['%{'];
            monLine{end+1}=['x_div = RAD2DEG*atan2(vx,vz);'];
            monLine{end+1}=['y_div = RAD2DEG*atan2(vy,vz);'];
            monLine{end+1}=['if (SCATTERED) flag=1; else flag=0;'];
            monLine{end+1}=['%}'];
            monLine{end+1}='';
            modifier=(i-1)*numel(monLine);
              
            
            B=[A(1:monitorPlacementList(i)+modifier) monLine A(monitorPlacementList(i)+modifier+1:end)];
            A=B;
        end
        
        % Insert in file
        fileID=fopen([folderPath '/' InstrumentName '/' InstrumentName '_analyze_stepwise.instr'],'w');
        for i=1:length(A)
            fprintf(fileID,'%s\n',char(A{i}));
        end
        fclose(fileID);
        fprintf('Created a stepwise analysis instrument with %i monitors along the guide \n',length(monitorPlacementList))
        % compile files
        if compileStepwise==0
            fileID=fopen([folderPath '/' InstrumentName '/compile_' InstrumentName '.sh'],'a');
            fprintf(fileID,'mcrun -n 0 -g -c --mpi ./%s_analyze_stepwise.instr\n',InstrumentName)
            fclose(fileID);
            fileID=fopen([folderPath '/' InstrumentName '/compile_' InstrumentName '_py.sh'],'a');
            fprintf(fileID,'mcrun-py -n 0 -g -c --mpi=1 %s_analyze_stepwise.instr\n',InstrumentName)
            fclose(fileID);
            compileStepwise=1;
        end
        % Change analyse file (CW)
end    

%% Add a optimize file without price if mode is not 'coating'
if strcmp(coatingOptions.scanType,'coating') == 0
    copyfile([coatingOptions.filePath slash InstrumentName '_optimize.instr'],[coatingOptions.filePath slash InstrumentName '_optimize_noPrice.instr'])
    fileID=fopen([coatingOptions.filePath slash InstrumentName '_optimize_noPrice.instr'],'r');
    tline = fgetl(fileID);
    clear A;
    i = 1;
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fileID);
        A{i} = tline;
    end
    A{end}= [];
    fclose(fileID)
    
    fileID=fopen([coatingOptions.filePath slash InstrumentName '_optimize_noPrice.instr'],'w');
    for i=1:length(A)
        if strfind(A{i},'Price=sectionPrices+TotalSubstratePrice)')
            tmp=[A{i}(1:strfind(A{i},'Price=sectionPrices+TotalSubstratePrice)')-1),'Price=1.00)']
            A{i}=tmp;
        end    
        fprintf(fileID,'%s\n',char(A{i}));
    end
    fclose(fileID);
    % compile files
    if compileNopirce==0
        fileID=fopen([folderPath '/' InstrumentName '/compile_' InstrumentName '.sh'],'a');
        fprintf(fileID,'mcrun -n 0 -g -c --mpi ./%s_optimize_noPrice.instr\n',InstrumentName)
        fclose(fileID);
        fileID=fopen([folderPath '/' InstrumentName '/compile_' InstrumentName '_py.sh'],'a');
        fprintf(fileID,'mcrun-py -n 0 -g -c --mpi=1 %s_optimize_noPrice.instr\n',InstrumentName)
        fclose(fileID);
        compileNopirce=1;
    end
    
end
%% Dublicate and create launch files if scan
if strcmp(coatingOptions.scanType,'pricescan')
    nFolders=length(coatingOptions.Price);
    thisDir=[cd];
    CreateNewPriceSim(thisDir,coatingOptions.Price); 
end

%% Add automatic view of log when using./launch_all
cd(thisFolder)
launchAllFile=textread('launch_all.sh','%s');
if strcmp(launchAllFile{end},'RunStatus.txt')
else
    fileID=fopen('launch_all.sh','a')
    fprintf(fileID,'watch -n 1 cat RunStatus.txt')
    fclose(fileID)
end




%% Save changes and print status
fprintf('\nDone inserting in %2.2f seconds. Files are ready to use\n',toc)

