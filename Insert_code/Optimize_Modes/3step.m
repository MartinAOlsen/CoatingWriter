options_single_home.ncount=1*1e7;
options_single_cluster.ncount=1*5e8;
options_single_cluster.ncount=1*1e8;


%% 3-step optimization process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This feature is added from CoatingWriter. See documentation for info %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Don't use this optimization process in "budget" or "price-scan" mode in
% CoatingWriter.
%
% This process runs three optimizations with different locked variables. 
% FIRST run:  Optimize only geometry, with locked coating distributins.
% SECOND run: Optimize only coating, with locked geometry.
% THIRD run:  Optimize only geometry with previous best pars as initial guess.
%
% If +A option is activated, an analysis will be performed after each
% optimization step. Be aware that this might take longer than 24 hours.
%
% If +C option is activated, the first run will be with m=6 and without 
% calculated value.
%
% Names of all the coatingvariables:
% INSERT COATING VARIABLES HERE!
pause(rand()*15)
printScript=fileread('printStatusCW.m');
tmpFolder=cd;
cd ..
outputFile=fileread('RunStatus.txt');
cd(tmpFolder)
runID=1;
while true
	if strfind(outputFile,sprintf('%3i |',runID))
		runID=runID+1
    else
        break
	end
end

method='3step';
tic;
steps=3;

%%%%%%%%%%
printStep='initialize';
time=toc;
eval(printScript)
%%%%%%%%%%

fileID=fopen('out_CW.txt','w')
fprintf(fileID,'Optimization Status by CoatingWriter\n')
fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf(fileID,'Optimizing using the %s method. \n',method)
fprintf(fileID,'Total steps in process: %i \n',steps)
fprintf(fileID,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')


fileID=fopen('out_CW.txt','w')
fclose(fileID)

ParNames=fieldnames(p);
try
    delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
end

% Make a list of every variable parameter for every run
CoatingParList={};
GeometryParList={};
AllParList={};
for i=1:length(ParNames) 
    if isa(eval(['p.' ParNames{i}]),'double')
        if any(strcmp(coatingVars,ParNames{i}))
            CoatingParList{end+1}=ParNames{i};
        else
            GeometryParList{end+1}=ParNames{i};
        end
        AllParList{end+1}=ParNames{i}; 
    end
end

%% First run (geometry):
Timing_FirstRun=tic;
step=1;
% lock variables to guess-values
tempPar=p;
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=num2str(p.' CoatingParList{i} '(2))'])
end

% INSERT FIRST INSTRUMENT HERE!
pars1=pars;
runResults(1)=o.criteriaBest;
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=ValueList;
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
filenameA=[filename '_geometry'];
scannameA=[scanname '_geometry'];
%% Update p guess
pFNames=fieldnames(pars);
for i = 1 : length(pFNames)
    eval(['p.' pFNames{i} '= [p.' pFNames{i} '(1), pars. ' pFNames{i} ',p.' pFNames{i} ']'])
end

TimeTable.firstRun=toc(Timing_FirstRun);

%% ADDANALYSIS %%


%% Second run (coating):
step=2;
Timing_SecondRun=tic;
tempPar=p;
for i=1:length(GeometryParList)
    eval(['tempPar.' GeometryParList{i} '=num2str(pars.' GeometryParList{i} ')'])
end

[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],tempPar,options{select});

runResults(2)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=[];
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
filenameA=[filename '_coating'];
scannameA=[scanname '_coating'];
%% Update p guess
pFNames=fieldnames(pars);
for i = 1 : length(pFNames)
    eval(['p.' pFNames{i} '= [p.' pFNames{i} '(1), pars. ' pFNames{i} ',p.' pFNames{i} ']'])
end
TimeTable.secondRun=toc(Timing_SecondRun);

%% ADDANALYSIS %%

%% Third run (geometry):
step=3;
tempPar=p;
Timing_ThirdRun=tic;
for i=1:length(GeometryParList)
    eval(['tempPar.' GeometryParList{i} '=[p.' GeometryParList{i} '(1) pars1.' GeometryParList{i} ' p.' GeometryParList{i} '(3)]'])
end

[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],tempPar,options{select});
runResults(3)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=[];
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
TimeTable.thirdRun=toc(Timing_ThirdRun);

filenameA=[filename '_geometry_2'];
scannameA=[scanname '_geometry_2'];
%% Update p guess
pFNames=fieldnames(pars);
for i = 1 : length(pFNames)
    eval(['p.' pFNames{i} '= [p.' pFNames{i} '(1), pars. ' pFNames{i} ',p.' pFNames{i} ']'])
end

%% ADDANALYSIS %%


Timing_StepScan=tic;
%% ADDSTEPSCAN %%
TimeTable.stepScan=toc(Timing_StepScan);


%%%%%%%%%% Print
printStep='analyze';
time=toc;
eval(printScript)
%%%%%%%%%%

monitor_ideal=monitor;
