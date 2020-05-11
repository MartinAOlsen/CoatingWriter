%% Begin StepScan
% This mode is a CoatingWriter feature. Look up latest documentation for
% further information.

% The function runs over each segment of the guide and tests every possible
% m-value. This works best after a CoatingWriter optimization using iFit,
% as a good initial guess is very important for this method to work.






%% Start status-plot if on desktop
supressPlots=0;
if select==1 && supressPlots==0
    try
	close all;
        FIG=figure('position',[1 1 1920 1080], 'color', 'w');
        subplot(2,3,3);
        title('Status')
        subplot(2,3,6);
        title('Last segment criteria')
        subplot(2,3,[1,2]);
        title('Coating distribution Horizontal')
        subplot(2,3,[4,5]);
        title('Coating distribution Vertical')
    catch
        supressPlots=1;
    end
end

%% print status
    totalSegments=0;
    doneSegments=0; % These values are used for the status UI
    printStep='stepScan';
    time=toc;
    try;eval(printScript);end

%% Set settings
    options_home.mode='simulate';
    options_cluster.mode='simulate';
    options_cluster.ncount=1*1e7;
    options_home.ncount=1*1e6;
    options_cluster.seed =1;
    options_home.seed =1;
    options_home = rmfield(options_home,'OutputFcn'); % Supress output from iFit
    options_home.Display='off';
    options={options_home options_cluster};
    

%% Lock parameters
    parlist=fieldnames(pars)
    for i=1:length(parlist)
        eval(['p.' parlist{i} '=''' num2str(eval(['pars.' parlist{i}])) ''';']);
    end
    % Set the stepScan variable to "1" changes the mode from normal CW/iFit to stepScan.
    p.stepScan='1';
    
%% Set up correction array
    % Load data from file
    % number of lines in file: 
    % get each line and put into a propper array:
    fid=fopen(['CoatingWriter_rawData' scanname '.txt']);
    n = 0;
    %tline = fgetl(fid);
    mValuesArray=[]; % Structure:  I , mTop , mBotom , mLeft , mRight
    tmpArray=[];
    lenArray=[];
    tmplenArray=[];
    while true
        tline = fgetl(fid);
        if ischar(tline)
            thisLine_splitted=textscan(tline,'%f,%1s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s');
            tmpArray(end+1,:)=[thisLine_splitted{1},thisLine_splitted{8},thisLine_splitted{9},thisLine_splitted{10},thisLine_splitted{11}];
            tmplenArray(end+1)=thisLine_splitted{3};
        else
            break
        end
    end
    fclose(fid)
    % Flip data 
    for i = max(tmpArray(:,1)):-1:1
        for j=1:length(tmpArray)
            if tmpArray(j,1) == i
                mValuesArray(end+1,:)=tmpArray(j,:);
                lenArray(end+1)=tmplenArray(j) ;
            end
        end
    end
    lenArray(1)=0;
    
    if select==1 && supressPlots==0
        accuLen=cumsum(lenArray)-0.5*lenArray;
        subplot(2,3,[1,2]);
        scatter(accuLen,mValuesArray(:,2),'k');
        axis([0,accuLen(end),0,6])
        hold on
        scatter(accuLen,mValuesArray(:,3),'k');
        axis([0,accuLen(end),0,6])
        hold off
        title('Coating distribution Horizontal (left and right mirror)')
        xlabel('Length from guide start [m]')
        ylabel('mvalue')
        
        
        subplot(2,3,[4,5]);
        scatter(accuLen,mValuesArray(:,4),'k');
        axis([0,accuLen(end),0,6])
        hold on
        scatter(accuLen,mValuesArray(:,5),'k');
        axis([0,accuLen(end),0,6])
        hold off
        title('Coating distribution Vertical (top and bottom mirror)')
        xlabel('Length from guide start [m]')
        ylabel('mvalue')
    end
    
    
    % see how many segments are in every section of the guide:
    segmentList=mValuesArray(:,1);
    counter=1;
    while true
        nOccurences(counter) = sum(segmentList==counter);
        if nOccurences(counter)==0 && counter > max(segmentList);
            nOccurences(counter)=[];
            break 
        end
        counter=counter+1;
    end
    
    %% Run once for to test
testRunOptions=options{select};
% testRunOptions.ncount=1e7;
[pars_test,monitor_test,m_test,o_test]=mcstas([instrument_name '_optimize.instr'],p,testRunOptions);
tmp=load(['priceListPunished' scanname '.txt']);
Price=tmp(end);
if select==1
    critBefore=monitor_test.Data.Criteria;
else
    critBefore=mean(mean(monitor(3).Data.I));  % Get the criteria like this on cluster, for some reason struct is different...
end
priceBefore=Price;
stepScanTime=tic;
    
    
%% Loop over sections
    
     for segment = length(nOccurences):-1:1  % optimize from moderator towards sample
        if isstr(fixSides{segment})
                switch fixSides{segment}
                case 'all'
                    IndexToRun{segment,1}={1,2,3,4};
                case 'none'
                    IndexToRun{segment,1}={1};
                    IndexToRun{segment,2}={2};
                    IndexToRun{segment,3}={3};
                    IndexToRun{segment,4}={4};
                case 'HV'
                    IndexToRun{segment,1}={1,2};
                    IndexToRun{segment,2}={3,4};
                case 'curve'
                    IndexToRun{segment,1}={1,2};
                    IndexToRun{segment,2}={3};
                    IndexToRun{segment,3}={4};
                end
        end
        totalSegments=totalSegments+size(IndexToRun,2)*nOccurences(segment);
    end
%% Loop over all subsegments (mirrors)
    for segment = length(nOccurences):-1:1  % Loop over each guide segment
        thisID = IndexToRun(segment,:);
        thisID(~cellfun('isempty',thisID))
        for index = 1 : sum(~cellfun('isempty', thisID( :)))
            mirrorsThisRun=IndexToRun{segment,index};
            
            for subSegment = 1 : nOccurences(segment) % Loop over each set of mirrors
                bestCrit=0;
                SegIndex=subSegment+sum(nOccurences(segment+1:end));
                ArrayIndex = sum(nOccurences)-SegIndex;
                %% Chose new m-value
                critList=[];
                errorList=[];
                for newMvalue = 1:0.5:5 % Loop over each possible m-value (TBD: make this varaible)
                    
                %% Update array
                    fid=fopen(['CoatingWriter_rawData' scanname '.txt'],'r');
                    n = 0;
%                     tline = fgetl(fid);
                    Array={};
                    mValuesArray=[]; % Structure:  I , mTop , mBotom , mLeft , mRight
                    while true
                        tline = fgetl(fid);
                        if ischar(tline)
                            thisLine_splitted=textscan(tline,'%f,%1s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s');
                            Array{end+1}=thisLine_splitted;
                        else
                            break
                        end
                    end
                    fclose(fid)
                    
                    % Insert new m-value:
                    tmp=Array{ArrayIndex};
                    for i=1:size(mirrorsThisRun,2)
                        tmp{mirrorsThisRun{i}+7}=newMvalue;
                    end
                    Array{ArrayIndex}=tmp;
                    
                    fid=fopen(['CoatingWriter_rawData' scanname '.txt'],'w');
                    for i=1:length(Array)
                       tmp=Array{i};
                       fprintf(fid,'%1.0f,%1s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n',tmp{1},char(tmp{2}),tmp{3},tmp{4},tmp{5},tmp{6},tmp{7},tmp{8},tmp{9},tmp{10},tmp{11},char(tmp{12}))
                    end
                    fclose(fid)
                    %% Run McStas
                    %[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],p,options{select});
                    [monitor] = simpleSim([instrument_name '_optimize.instr'],p,options{select});
                    iters=iters+1;
                
                    %% Save part-result
                    % Calculate criteria:
                    if select==1
%                         Criteria=monitor.Data.Criteria;
%                         Error=monitor.Data.values(2);
                            Criteria=monitor.Data.Criteria
                            Error = monitor.Data.Error
                    else
%                         Criteria=mean(mean(monitor(3).Data.I)); 
%                         Error=monitor.Data.values(2);
                            Criteria=monitor.Data.Criteria
                            Error = monitor.Data.Error
                    end
                    % If criteria is new best, save
                    if abs(Criteria)>abs(bestCrit)
                        bestCorr=newMvalue;
                        bestCrit=Criteria;
                        tmp=load(['priceListPunished' scanname '.txt']);
                        Price=tmp(end);
                        
                    end
                    delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
                    critList(end+1)=Criteria;
                    errorList(end+1)=Error;
                    
                    if select==1 && supressPlots==0
                        %figure(FIG)
                        % Plot last section crits
                        subplot(2,3,6);
                        scatter(1:0.5:newMvalue,abs(critList),'b','filled');
                        hold on
                        %errorbar(1:0.5:newMvalue,abs(critList),abs(errorList),'LineStyle','none','Color','r')
			pause(0.05)

                        xlabel('m-value')
                        ylabel('criteria (higher is better)')
                        title('Current scan')
                    end
                
                end
                
                
                iters=iters+1;
                doneSegments=doneSegments+1;
                %% Update status
                 time=toc;
                 try;eval(printScript);end
                %% Update array (after scan)
                    fid=fopen(['CoatingWriter_rawData' scanname '.txt'],'r');
                    n = 0;
%                     tline = fgetl(fid);
                    Array={};
                    mValuesArray=[]; % Structure:  I , mTop , mBotom , mLeft , mRight
                    while true
                        tline = fgetl(fid);
                        if ischar(tline)
                            thisLine_splitted=textscan(tline,'%f,%1s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s');
                            Array{end+1}=thisLine_splitted;
                            mValuesArray(end+1,:)=[thisLine_splitted{1},thisLine_splitted{8},thisLine_splitted{9},thisLine_splitted{10},thisLine_splitted{11}];
                        else
                            break
                        end
                    end
                    fclose(fid)
                    
                    % Insert new m-value:
                    tmp=Array{ArrayIndex};
                    for i=1:size(mirrorsThisRun,2)
                        tmp{mirrorsThisRun{i}+7}=bestCorr;
                    end
                    Array{ArrayIndex}=tmp;
                    
                    fid=fopen(['CoatingWriter_rawData' scanname '.txt'],'w');
                    for i=1:length(Array)
                       tmp=Array{i};
                       mValuesArray(i,:)=[tmp{1},tmp{8},tmp{9},tmp{10},tmp{11}];
                       fprintf(fid,'%1.0f,%1s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n',tmp{1},char(tmp{2}),tmp{3},tmp{4},tmp{5},tmp{6},tmp{7},tmp{8},tmp{9},tmp{10},tmp{11},char(tmp{12}))
                    end
                    fclose(fid) 

            %% Try to plot progress
                    if select==1 && supressPlots==0
                        %figure(FIG)
                    % Plot last section crits
                        subplot(2,3,6);
                        scatter(1:0.5:5,abs(critList),'b','filled');
                        hold on
                        %errorbar(1:0.5:5,abs(critList),abs(errorList),'LineStyle','none','Color','r')
                        hold off
                        xlabel('m-value')
                        ylabel('criteria (higher is better)')
                        title('Current scan')
                        
                    % Update table of values
                        hl1 = subplot(2,3,3);
                        cla(hl1)
                        text(0.1,1,'CoatingWriter StepScan','FontSize',22)
                        text(0,0.83,sprintf('%16s: %3.1f%%','Progress',100*doneSegments/totalSegments),'FontSize',16,'HorizontalAlignment','left')
                        if Price-priceBefore > 0; isPos='+'; else; isPos=''; end
                        text(0,0.66,sprintf('%16s: %4.1f (%s%3.1f%%)','Price',Price,isPos,100*(Price-priceBefore)/priceBefore),'FontSize',16,'HorizontalAlignment','left');
                        if strcmp(Result_this.mode,'coating')
                             if bestCrit-critBefore > 0; isPos='+'; else; isPos=''; end
                             text(0,0.5,sprintf('%16s: %4.3f (%s%3.1f%%)',Intensity,bestCrit,isPos,100*(bestCrit-critBefore)/critBefore),'FontSize',16);
                        else
                            I_now=abs(bestCrit*Price);
                            I_before=abs(critBefore*priceBefore);
                            if I_now-I_before > 0; isPos='+'; else; isPos=''; end
                            text(0,0.5,sprintf('%16s: %4.3f (%s%3.1f%%)','Intensity',I_now,isPos,100*(I_now-I_before)/I_before),'FontSize',16,'HorizontalAlignment','left');
                        end
                        text(0,0.33,sprintf('%16s: %2.1f minutes','Simulation time',toc(stepScanTime)/60),'FontSize',16,'HorizontalAlignment','left')
                        
                        timePrPercent=toc(stepScanTime)/(100*doneSegments/totalSegments);
                        timeLeft=timePrPercent*(100-100*doneSegments/totalSegments);
                        text(0,0.17,sprintf('%16s: %2.1f minutes','Time left',timeLeft/60),'FontSize',16,'HorizontalAlignment','left')
                        
                        
                        ax = gca
                        ax.Visible = 'off';
%                         uitable('Data',{1,2,'a'})
                    
                    % Plot New distribution
                        %thisSegmentNumber=mValuesArray(ArrayIndex,1);
                        %accuLen=cumsum(lenArray);
                        %nOffset=sum(nOccurences(end:-1:segment+1));
                        %arrayOffset=sum(nOccurences(1:segment-1));
                        %lenOffset=sum(lenArray(1:nOffset));
                        %thisSegmentOffset=cumsum(lenArray(nOffset+1:ArrayIndex-arrayOffset+nOffset));
                        %if isempty(thisSegmentOffset) ; thisSegmentOffset = 0; end
                        %thisLength=thisSegmentOffset(end)+lenOffset-0.5*(lenArray(ArrayIndex-arrayOffset));
                         thisLength= sum(lenArray(1:SegIndex-1)) + 0.5* lenArray(SegIndex);
                        
                         mirrorToPlot=[mirrorsThisRun{:}];
                        
                        subplot(2,3,[1,2]);
                        hold on
                        if sum(mirrorToPlot==1)>0
                            %scatter(thisLength,mValuesArray(ArrayIndex,2)','g','filled');
                            line([thisLength-0.5* lenArray(SegIndex),thisLength+0.5* lenArray(SegIndex)],[mValuesArray(ArrayIndex,2),mValuesArray(ArrayIndex,2)],'color','green','linewidth',3)
                        end
                        
                        axis([0,accuLen(end),0,6])
                        if sum(mirrorToPlot==2)>0
                            %scatter(thisLength,mValuesArray(ArrayIndex,3)','b','filled');
                            line([thisLength-0.5* lenArray(SegIndex),thisLength+0.5* lenArray(SegIndex)],[mValuesArray(ArrayIndex,3),mValuesArray(ArrayIndex,3)],'color','blue','linewidth',3)
                        end
                        hold off
                        axis([0,accuLen(end),0,6])
                        subplot(2,3,[4,5]);
                        hold on
                        if sum(mirrorToPlot==3)>0
                            line([thisLength-0.5* lenArray(SegIndex),thisLength+0.5* lenArray(SegIndex)],[mValuesArray(ArrayIndex,4),mValuesArray(ArrayIndex,4)],'color','magenta','linewidth',3)
                            %scatter(thisLength,mValuesArray(ArrayIndex,4)','m','filled');
                        end
                        
                        axis([0,accuLen(end),0,6])
                        if sum(mirrorToPlot==4)>0
                            line([thisLength-0.5* lenArray(SegIndex),thisLength+0.5* lenArray(SegIndex)],[mValuesArray(ArrayIndex,5),mValuesArray(ArrayIndex,5)],'color','red','linewidth',3)
                            %scatter(thisLength,mValuesArray(ArrayIndex,5)','r','filled');
                        end
                        axis([0,accuLen(end),0,6])
                        hold off
                        title('Coating distribution Vertical')
                        pause(0.15)
                    end
            
            end
        end
    end
    TotalStepScanTime=toc(stepScanTime);
    
%% Save result
