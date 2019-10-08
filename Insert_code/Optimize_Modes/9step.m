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
% FIRST run:   Optimize geometry.
% SECOND run:  Optimize coating with previous best pars as initial guess.
% THIRD run:   Optimize all with previous best pars as initial guess.
% FOURTH run:  Optimize geometry with previous best pars as initial guess
% FIFTH run:   Optimize coating with previous best pars as initial guess.
% SIXTH run:   Optimize all with previous best pars as initial guess.
% SEVENTH run: Optimize geometry with previous best pars as initial guess
% EIGHTH run:  Optimize coating with previous best pars as initial guess.
% NINETH run:  Optimize all with previous best pars as initial guess.
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
method='6step';
steps=9;
tic
%%%%%%%%%%
printStep='initialize';
time=toc;
eval(printScript)
%%%%%%%%%%

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
% lock variables to guess-values
step=1;
tempPar=p;
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=num2str(p.' CoatingParList{i} '(2))'])
end

% INSERT FIRST INSTRUMENT HERE!
pars1=pars;
runResults(1)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:end-1);
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=ValueList;
%%%%%%%%%%
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);


filenameA=[filename '_geometry'];
scannameA=[scanname '_geometry'];


%% ADDANALYSIS %%


%% Second run (coating):
step=2;
tempPar=p;
for i=1:length(GeometryParList)
    eval(['tempPar.' GeometryParList{i} '=num2str(pars.' GeometryParList{i} ')'])
end

[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],tempPar,options{select});
runResults(2)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=[];
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:end-1);
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
%%%%%%%%%%
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
filenameA=[filename '_coating'];
scannameA=[scanname '_coating'];


%% ADDANALYSIS %%

%% Third run (all):
step=3;
tempPar=p;
for i=1:length(GeometryParList)
    eval(['tempPar.' GeometryParList{i} '=[p.' GeometryParList{i} '(1) pars1.' GeometryParList{i} ' p.' GeometryParList{i} '(3)]'])
end
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=[p.' CoatingParList{i} '(1) pars.' CoatingParList{i} ' p.' CoatingParList{i} '(3)]'])
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
%printStatusCW(o,ValueList,mode,3,toc,'6step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
filenameA=[filename '_all'];
scannameA=[scanname '_all'];
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%

for i=1:length(GeometryParList)
    eval(['p.' GeometryParList{i} '=[p.' GeometryParList{i} '(1) pars.' GeometryParList{i} ' p.' GeometryParList{i} '(3)]'])
end
for i=1:length(CoatingParList)
    eval(['p.' CoatingParList{i} '=[p.' CoatingParList{i} '(1) pars.' CoatingParList{i} ' p.' CoatingParList{i} '(3)]'])
end

%% ADDANALYSIS %%


%% Fourth run (geometry):
step=4;
% lock variables to guess-values
tempPar=p;
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=num2str(p.' CoatingParList{i} '(2))'])
end

[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],tempPar,options{select});

pars1=pars;
runResults(1)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
%%%printStatusCW(o,ValueList,mode,4,toc,'9step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);


filenameA=[filename '_geometry_2'];
scannameA=[scanname '_geometry_2'];
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
%% ADDANALYSIS %%


%% Fifth run (coating):
step=5;
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
%printStatusCW(o,ValueList,mode,5,toc,'9step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
filenameA=[filename '_coating_2'];
scannameA=[scanname '_coating_2'];
% Find last values (crit,I,price)
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
%% ADDANALYSIS %%

%% Sixth run (all):
tempPar=p;
step=6;
for i=1:length(GeometryParList)
    eval(['tempPar.' GeometryParList{i} '=[p.' GeometryParList{i} '(1) pars1.' GeometryParList{i} ' p.' GeometryParList{i} '(3)]'])
end
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=[p.' CoatingParList{i} '(1) pars.' CoatingParList{i} ' p.' CoatingParList{i} '(3)]'])
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
%printStatusCW(o,ValueList,mode,6,toc,'9step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
filenameA=[filename '_all_2'];
scannameA=[scanname '_all_2'];
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
for i=1:length(GeometryParList)
    eval(['p.' GeometryParList{i} '=[p.' GeometryParList{i} '(1) pars.' GeometryParList{i} ' p.' GeometryParList{i} '(3)]'])
end
for i=1:length(CoatingParList)
    eval(['p.' CoatingParList{i} '=[p.' CoatingParList{i} '(1) pars.' CoatingParList{i} ' p.' CoatingParList{i} '(3)]'])
end
%% ADDANALYSIS %%



%% Seventh run (geometry):
step=7;
% lock variables to guess-values
tempPar=p;
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=num2str(p.' CoatingParList{i} '(2))'])
end

[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],tempPar,options{select});

pars1=pars;
runResults(1)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=load(['priceList' scanname '.txt']);
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load(['priceListPunished' scanname '.txt']);
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
%printStatusCW(o,ValueList,mode,7,toc,'9step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);


filenameA=[filename '_geometry_3'];
scannameA=[scanname '_geometry_3'];
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
%% ADDANALYSIS %%


%% Eighth run (coating):
tempPar=p;
step=8;
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
%printStatusCW(o,ValueList,mode,8,toc,'9step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
filenameA=[filename '_coating_3'];
scannameA=[scanname '_coating_3'];
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
%% ADDANALYSIS %%

%% Nineth run (all):
step=9;
tempPar=p;
for i=1:length(GeometryParList)
    eval(['tempPar.' GeometryParList{i} '=[p.' GeometryParList{i} '(1) pars1.' GeometryParList{i} ' p.' GeometryParList{i} '(3)]'])
end
for i=1:length(CoatingParList)
    eval(['tempPar.' CoatingParList{i} '=[p.' CoatingParList{i} '(1) pars.' CoatingParList{i} ' p.' CoatingParList{i} '(3)]'])
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
%printStatusCW(o,ValueList,mode,9,toc,'9step')
delete (['priceList' scanname '.txt']);delete (['priceListPunished' scanname '.txt']);
filenameA=[filename '_all_3'];
scannameA=[scanname '_all_3'];
%%%%%%%%%% Print
printStep='optimize';
time=toc;
eval(printScript)
%%%%%%%%%%
%% ADDANALYSIS %%


%% ADDSTEPSCAN %%


%%%%%%%%%% Print
printStep='analyze';
time=toc;
eval(printScript)
%%%%%%%%%%
monitor_ideal=monitor;
