%% 3-step optimization process with single parameter scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This feature is added from CoatingWriter. See documentation for info %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Don't use this optimization process in "budget" or "price-scan" mode in
% CoatingWriter.
%
% This process runs three optimizations with different locked variables. 
% FIRST run:  Optimize all variables.
% SECOND run: Optimize all variables one at a time.
% THIRD run:  Optimize all variables previous best pars as initial guess.
Result_this.OptimizeMode='singleParScan';


tic
fileID=fopen('out_CW.txt','w')
fclose(fileID)

ParNames=fieldnames(p);
try
    delete priceList.txt;delete priceListPunished.txt;
end
ParList=[];
for i=1:length(ParNames) 
    if isa(eval(['p.' ParNames{i}]),'double')
	ParList{end+1}=ParNames{i};
    end
end

%% First run:
step=1;
[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],p,options{select});
pars1=pars;
runResults(1)=o.criteriaBest;
Result_this.critList(step)=o.criteriaBest;
ValueList=load('priceList.txt');
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load('priceListPunished.txt');
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=ValueList;
%printStatusCW(o,ValueList,mode,1,toc,'singleParScan')
fileID=fopen('out_CW.txt','a')
fprintf(fileID,'Step 1/3\n\tAll parameters optimized\n\t\tCrit = %d\n',o.criteriaBest)
fclose(fileID)

delete priceList.txt;delete priceListPunished.txt;
for i=1:length(ParList)
    eval(['p.' ParList{i} '=[p.' ParList{i} '(1),pars.' ParList{i} ',p.' ParList{i} '(3)]'])
end
%% ADDANALYSIS %%

%% Second run single parameter scan:
step=2;
tempOptions=options{select};
tempOptions.MaxFunEvals=100;
tempOptions.MaxIter=100;
tempOptions.TolFun ='0.15%';
tempOptions.TolX ='0.15%';
tempOptions.optimizer='fminimga';
History=[];
fileID=fopen('out_CW.txt','a')
fprintf(fileID,'Step 2/3\n\tSingle Par scan:\n')
fclose(fileID)
oldCrit=o.criteriaBest;

for i=1:length(ParList)
	% set parameters
	tempPar=p;
	for j=1:length(ParList)
	    if i==j
    	        
	    else
		eval(['tempPar.' ParList{j} '=num2str(p.' ParList{j} '(2));'])
	    end
	end
	[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],tempPar,tempOptions);
	eval(['p.' ParList{i} '=[p.' ParList{i} '(1), o.parsBest ,p.' ParList{i} '(3)];'])
	History=[History;o.criteriaHistory];
	runResults=o.criteriaBest;
	Result_this.critList(end+1)=runResults;
	fileID=fopen('out_CW.txt','a')
	fprintf(fileID,'\t\t Optimized %s, new crit: %d, changed %2.2f%%\n',ParList{i},o.criteriaBest,100*(o.criteriaBest-oldCrit)/oldCrit)
	fclose(fileID)
	oldCrit=o.criteriaBest;
end


ValueList=[];
ValueList=load('priceList.txt');
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load('priceListPunished.txt');
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
%printStatusCW(o,ValueList,mode,2,toc,'singleParScan')
fileID=fopen('out_CW.txt','a')
fprintf(fileID,'Step 3/3\n\tAll parameters optimized\n\t\tCrit = %d\n',o.criteriaBest)
fclose(fileID)
delete priceList.txt;delete priceListPunished.txt;
%% ADDANALYSIS %%

%% Third run (all):
step=3;
[pars,monitor,m,o]=mcstas([instrument_name '_optimize.instr'],p,options{select});
runResults=o.criteriaBest;
Result_this.critList(end+1)=runResults;
ValueList=[];
ValueList=load('priceList.txt');
ValueList=ValueList(1:length(o.criteriaHistory));
ValueList(:,2)=o.criteriaHistory';
tmp=load('priceListPunished.txt');
ValueList(:,3)=tmp(1:length(o.criteriaHistory));
Result_this.valueList=[Result_this.valueList;ValueList];
%printStatusCW(o,ValueList,mode,3,toc,'singleParScan')
delete priceList.txt;delete priceListPunished.txt;
for i=1:length(ParList)
    eval(['p.' ParList{i} '=[p.' ParList{i} '(1),pars.' ParList{i} ',p.' ParList{i} '(3)]'])
end
%% ADDANALYSIS %%
