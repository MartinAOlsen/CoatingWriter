function [Lines_ifit,Lines_instr]=CoatingLineWriter_Value(coatingOptions,i)
%COATINGLINEWRITER_VALUE2 Summary of this function goes here
%   second build of this function, first introduced in CoatingWriter 2.1.05
%   created on 26/04-2017 by Martin Olsen

% Takes input from CoatingWriter, and returns the lines for guide_bot-files


numSecPrElement=500;
if isfield(coatingOptions,'mStepLength') && isfield(coatingOptions,'instrumentLength')
    numSecPrElement=round(coatingOptions.instrumentLength/coatingOptions.mStepLength)+1;
end


%% Initialize arrays
tic
Print_Bisection=0;
Lines_ifit={};
listOfInput=fieldnames(coatingOptions);
C_print=1;
Lines_instr.COATINGfirst={};
Lines_instr.COATING={};
Lines_instr.DEFINE={};
Lines_instr.DECLARE={};
Lines_instr.Monitor={};
Lines_instr.SUMMARY={};
Lines_instr.SUMMARY={};
Lines_instr.INITIALIZE={};
Lines_instr.DEFINE{1}='//Define'
Lines_instr.DECLARE{1}='//Declare'
Declare={};
Define={};
Defined={};
Declared={};

%% Some declerations
if i==1; Lines_instr.DECLARE{end+1}=sprintf('double localCounter;');end
if i==1; Lines_instr.DECLARE{end+1}=sprintf('double Price;');end
Declare{end+1}=sprintf('Segment%iSubstrateSize',i)
Define{end+1}=sprintf('beta')
Declare{end+1}=sprintf('sumOfMValues')
Declare{end+1}=sprintf('sectionPrices')
Declare{end+1}=sprintf('power')
if i==1;Declare{end+1}=sprintf('priceToSave');end
if i==1
   Declare{end+1}=sprintf('verticalSegmentSize [%i]',numSecPrElement);
   Declare{end+1}=sprintf('horizontalSegmentSize [%i]',numSecPrElement);     
end

Define{end+1}=sprintf('guidePrice')

%% M-step options 
mStepOption=[sum(strcmp(listOfInput,'deltaM')),sum(strcmp(listOfInput,'mStepLength'))];
if mStepOption(1)>1 || mStepOption(2)>1
    error('Something is wrong with the m-value steps. More than one either deltaM or mStepLength was found in input. Please address this!')
end
if mStepOption(1)==1 && mStepOption(2)==1
    useMSteps=1;
else
    useMSteps=0;
end

%% Begin
fprintf('Making coating parameters for segment %i (%s) with CoatingLineWriter_Value.m',i,coatingOptions.segmentType)

mStepOption=[sum(strcmp(listOfInput,'deltaM')),sum(strcmp(listOfInput,'mStepLength'))];

if i==1 
   Lines_instr.DEFINE{end+1}='stepScan=0,';  %% Should be =0
end
Lines_ifit{end+1}=['fixSides{' num2str(i) '} = ''' coatingOptions.fixSides ''''];



%% Determine number of free panels
if strcmp(coatingOptions.fixSides,'all')
    freePanels=1;
elseif strcmp(coatingOptions.fixSides,'horizontal') || strcmp(coatingOptions.fixSides,'Horizontal') || strcmp(coatingOptions.fixSides,'horisontal')
    freePanels=2;
    freePanelList=[34,1,2];
elseif strcmp(coatingOptions.fixSides,'vertical')
    freePanels=2;
    freePanelList=[12,3,4];
elseif strcmp(coatingOptions.fixSides,'none')
    freePanels=4;
    freePanelList=[1,2,3,4];
elseif strcmp(coatingOptions.fixSides,'curve')
    freePanels=3;
    freePanelList=[12,5,6];
elseif strcmp(coatingOptions.fixSides,'HV') || strcmp(coatingOptions.fixSides,'VH')
    freePanels=2;
    freePanelList=[34,12];
end

if i==1
        Lines_ifit{end+1}=sprintf('if length(scanname)>0'); 
	Lines_ifit{end+1}=sprintf('\tp.scanname=scanname;'); 
	Lines_ifit{end+1}=sprintf('else'); 
	Lines_ifit{end+1}=sprintf('\tp.scanname=''noScanName'';'); 
	Lines_ifit{end+1}=sprintf('end'); 
    Lines_instr.DEFINE{end+1}=sprintf('string scanname="",');
end
%% Loop for every panel
for n=1:freePanels
    if freePanels==1; panel='';
    elseif (freePanelList(n)==1); panel='top';
    elseif (freePanelList(n)==2); panel='bottom';
    elseif (freePanelList(n)==3); panel='left';
    elseif (freePanelList(n)==4); panel='right';
    elseif (freePanelList(n)==5); 
        if strcmp(coatingOptions.type,'C')
            panel='inside';
        else
            panel='left';
        end
    elseif (freePanelList(n)==6); 
        if strcmp(coatingOptions.type,'C')
            panel='outside';
        else
            panel='right';
        end
    elseif (freePanelList(n)==34); panel='horizontal';
    elseif (freePanelList(n)==12); panel='vertical'; 
    end
    
    
    minA=coatingOptions.minA;
    maxA=coatingOptions.maxA;
    minB=coatingOptions.minB;
    maxB=coatingOptions.maxB;
    minC=coatingOptions.minC;
    maxC=coatingOptions.maxC;
    if strcmp(coatingOptions.distribution,'constant')
        coatingOptions.formula='';
    end
    this_formula=coatingOptions.formula;
    minCenter=coatingOptions.minCenter;
    maxCenter=coatingOptions.maxCenter;
    
    % Insert variable numbers.
    C_number=sprintf('c%d%s',i,panel);
    this_formula=strrep(this_formula,'C',C_number);
    poly_number=sprintf('Power%d%s',i,panel);
    this_formula=strrep(this_formula,'power',poly_number);
    Define{end+1}=sprintf('Power%d%s',i,panel);
    A_number=sprintf('a%d%s',i,panel);
    this_formula=strrep(this_formula,'A',A_number);
    B_number=sprintf('b%d%s',i,panel);
    this_formula=strrep(this_formula,'B',B_number);
    center_number=sprintf('center%d%s',i,panel);
    this_formula=strrep(this_formula,'center',center_number);
    
    % Number of segments. Make variable?
    segments=coatingOptions.nSegments;
    segments=500;
    
    
    if i==1&&n==1;Lines_ifit{end+1}=sprintf('%% CoatingWriter Options below:');else;end
%     if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.beta=''2.6''; %%Price exponent - Prices use: Price=(m*length)^Beta.');else;end
%     if i==1&&n==1;Lines_ifit{end+1}=sprintf('p.guidePrice=''0''; %%Base price pr meter of guide without coating.');else;end
    
    if n==1;
    if strcmp(coatingOptions.fixSides,'vertical')
        freePanelList=[12,3,4];
    elseif strcmp(coatingOptions.fixSides,'horizontal')
        freePanelList=[34,1,2];
    elseif strcmp(coatingOptions.fixSides,'none')
        freePanelList=[1,2,3,4];
    end
    end
    if freePanels==1
        numberInThisRun=1;
   else
         if (freePanelList(n)==1); numberInThisRun=1;
         elseif (freePanelList(n)==2); numberInThisRun=1;
         elseif (freePanelList(n)==3); numberInThisRun=1;
         elseif (freePanelList(n)==4); numberInThisRun=1;
         elseif (freePanelList(n)==5); numberInThisRun=1;
         elseif (freePanelList(n)==6); numberInThisRun=1;
         elseif (freePanelList(n)==34); numberInThisRun=2;
         elseif (freePanelList(n)==12); numberInThisRun=2; 
         elseif (freePanelList(n)==1234); numberInThisRun=4; 
         end
    end
    
    if strcmp(coatingOptions.distribution,'power') || strcmp(coatingOptions.distribution,'power')
        inputvariables={'a','b','Power','center'}
    elseif strcmp(coatingOptions.distribution,'linear')
        inputvariables={'a','b','center'}
    elseif strcmp(coatingOptions.distribution,'constant')
        inputvariables={'mValues'}
    elseif strcmp(coatingOptions.distribution,'exp')
        inputvariables={'a','b','c','center'}
    end
%     Declare{end+1}=sprintf('power%d%s',i,panel);
    
    
    %% Determine parameter limits
    for j=1:length(inputvariables)
        tmpVariable=char(inputvariables{j});
        if length(tmpVariable)==1
            tmpVariable=upper(tmpVariable);
        else
            if strcmp(tmpVariable,'Power')==1
                tmpVariable='Power';
            elseif strcmp(tmpVariable,'mValues')==1
                tmpVariable='M';
            else
                tmpVariable=[upper(tmpVariable(1)) tmpVariable(2:end)];
            end
        end

        tmpMin=eval(['coatingOptions.min' tmpVariable]);
        tmpMax=eval(['coatingOptions.max' tmpVariable]);

        if tmpMin==tmpMax
            Lines_ifit{end+1}=sprintf('p.%s%i%s=''%2.2f'';',char(inputvariables{j}),i,panel,tmpMax);
        else
            if strcmp(coatingOptions.guessType,'m6')
                localGuess=eval(sprintf('coatingOptions.guess.%s',char(inputvariables{j})));
                Lines_ifit{end+1}=sprintf('p.%s%i%s=[%2.2f %2.2f %2.2f];',char(inputvariables{j}),i,panel,tmpMin,localGuess,tmpMax); 
            elseif strcmp(coatingOptions.guessType,'mean')
                Lines_ifit{end+1}=sprintf('p.%s%i%s=[%2.2f %2.2f %2.2f];',char(inputvariables{j}),i,panel,tmpMin,(tmpMax+tmpMin)/2,tmpMax); 
            end
        end
        Define{end+1}=sprintf('%s%i%s',char(inputvariables{j}),i,panel);
    end
    
    
    %% Power distribution
    if strcmp(coatingOptions.distribution,'power') || strcmp(coatingOptions.distribution,'power')
        

        
        %%%%% START FOR LOOP
        Declare{end+1}=sprintf('Power')
        Declare{end+1}='counter2';
        Declare{end+1}=sprintf('mValues%i%s',i,panel)
        Lines_instr.INITIALIZE{1}=sprintf('////Segment %d',i);
        if useMSteps==1
            Lines_instr.DECLARE{end+1}=sprintf('int segment;')
            Lines_instr.COATING{end+1}=sprintf('for (segment =1 ; segment < nSegments%i+1 ; segment++)',i);
            Lines_instr.COATING{end+1}=sprintf('{');
            x=sprintf('(%s/Nsegments%i)-','localCounter',i);
        else
            Lines_instr.COATING{end+1}=sprintf('for (segment =1 ; segment < %i ; segment++)',segments+1);
            Lines_instr.COATING{end+1}=sprintf('{');
            x=sprintf('(%s/%d)-','localCounter',segments);
        end
        local_formula=strrep(this_formula,'x-',x);

        Lines_instr.COATING{end+1}=sprintf('\tmValues%d%s[segment-1] = %s ;',i,panel,local_formula);
        Lines_instr.DECLARE{end+1}=sprintf('double mValues%d%s [%i];',i,panel,numSecPrElement);
        Declared{end+1}=sprintf('mValues%d%s',i,panel);
        Lines_instr.DECLARE{end+1}=sprintf('double elementLength%d [%i];',i,numSecPrElement);
        Declared{end+1}=sprintf('elementLength%d',i);
%         tmpM=strrep(coatingOptions.alpha,'m',sprintf('mValues%d%s[segment-1] ',i,panel))
%         Lines_instr.COATING{end+1}=sprintf('\talpha%d%s[segment-1] = %s ;',i,panel,tmpM);
%         Lines_instr.DECLARE{end+1}=sprintf('double alpha%d%s [%i];',i,panel,segments);
        if useMSteps==1
            Lines_instr.COATING{end+1}=sprintf('    mValues%i%s[segment-1] = round(mValues%i%s[segment-1]*%2.4f)/%2.4f ;',i,panel,i,panel,1/coatingOptions.deltaM,1/coatingOptions.deltaM);
            Lines_instr.COATING{end+1}=sprintf('    elementLength%d[segment-1] = length%d/Nsegments%i ;',i,i,i);
        else
            Lines_instr.COATING{end+1}=sprintf('    elementLength%d[segment-1] = length%d/%d ;',i,i,segments);
        end
        if coatingOptions.remove_m_under > 0
            Lines_instr.COATING{end+1}=sprintf('\tif (mValues%i%s[segment-1] < %d)',i,panel,coatingOptions.remove_m_under);
            Lines_instr.COATING{end+1}=sprintf('\t\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\tmValues%i%s[segment-1] = %d;',i,panel,coatingOptions.remove_m_under);
            Lines_instr.COATING{end+1}=sprintf('\t\t}');
        end

        if coatingOptions.remove_m_over > 0
            Lines_instr.COATING{end+1}=sprintf('\tif (mValues%i%s[segment-1] > %d)',i,panel,coatingOptions.remove_m_over);
            Lines_instr.COATING{end+1}=sprintf('\t\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\tmValues%i%s[segment-1] = %d;',i,panel,coatingOptions.remove_m_over);
            Lines_instr.COATING{end+1}=sprintf('\t\t}');
        end
            % stepScan overwrite
            Lines_instr.COATING{end+1}=sprintf('\tif (stepScan == 1)');
            Lines_instr.COATING{end+1}=sprintf('\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\tif (correctionMValue%i%s[segment-1] > 0.1)',i,panel);
            Lines_instr.COATING{end+1}=sprintf('\t{\t');
            Lines_instr.COATING{end+1}=sprintf('\t\t\tmValues%i%s[segment-1] = correctionMValue%i%s[segment-1];',i,panel,i,panel);
            Lines_instr.DECLARE{end+1}=sprintf('double correctionMValue%i%s [1000];',i,panel);
            Lines_instr.COATING{end+1}=sprintf('\t}\t');
            Lines_instr.COATING{end+1}=sprintf('\t}');


            

%         if freePanels > 1
%             Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (elementLength%d[segment-1] * (pow(mValues%i%s[segment-1],beta))/%i)*%i ;',i,i,panel,4,numberInThisRun);  %% Sum of m-value-meter
%         else
%             Lines_instr.COATING{end+1}=sprintf('    sumOfMValues += (elementLength%d[segment-1] * pow(mValues%i%s[segment-1],beta)) ;',i,i,panel);  %% Sum of m-value-meter
%         end
        % Coating Price
        [Lines_instr,Lines_ifit] = coatingPrice(i,n,panel,coatingOptions,Lines_ifit,Lines_instr,numberInThisRun,freePanels)
        
        
        Lines_instr.COATING{end+1}=sprintf('localCounter++;');
        Lines_instr.COATING{end+1}=sprintf('}');
        Lines_instr.COATING{end+1}=sprintf('localCounter=0;');

        Lines_instr.COATING{end+1}=sprintf('');

%         Lines_instr.COATING{end+1}=sprintf('sectionPrices+=sumOfMValues;');
%         Lines_instr.COATING{end+1}=sprintf('sumOfMValues=0;');
        
        if useMSteps==1
		Lines_instr.COATING{end+1}='MPI_MASTER('
                Lines_instr.COATING{end+1}=sprintf('printf("\\n-----m-value %s (Start)=%%2.3f ----",mValues%i%s[0] );',panel,i,panel);
                %Lines_instr.COATING{end+1}=sprintf('printf("\\n-----m-value %s (Mid)=%%2.3f ----",mValues%i%s[round(Nsegments%i/2)] );',panel,i,panel,i)
                Lines_instr.COATING{end+1}=sprintf('printf("\\n-----m-value %s (End)=%%2.3f ----",mValues%i%s[nsegments%i-1] );',panel,i,panel,i);
		Lines_instr.COATING{end+1}=')'
        else 
		Lines_instr.COATING{end+1}='MPI_MASTER('
                Lines_instr.COATING{end+1}=sprintf('printf("\\n-----m-value %s (Start)=%%2.3f ----",mValues%i%s[0] );',panel,i,panel);
                Lines_instr.COATING{end+1}=sprintf('printf("\\n-----m-value %s (Mid)=%%2.3f ----",mValues%i%s[49] );',panel,i,panel);
                Lines_instr.COATING{end+1}=sprintf('printf("\\n-----m-value %s (End)=%%2.3f ----",mValues%i%s[99] );',panel,i,panel);
		Lines_instr.COATING{end+1}=')'
        end
        
        
        
        
        
    elseif strcmp(coatingOptions.distribution,'constant') 
%         Lines_instr.COATING{end+1}=sprintf('sumOfMValues += (length%i * (pow(mValues%i%s,beta))/%i)*%i ;',i,i,panel,4,numberInThisRun);  %% Sum of m-value-meter
%         Lines_instr.COATING{end+1}=sprintf('sectionPrices+=sumOfMValues;');
%         Lines_instr.COATING{end+1}=sprintf('sumOfMValues=0;');
%         Declare{end+1}=sprintf('mValues%i%s',i,panel)
        Declare{end+1}=sprintf('mValues%i%s',i,panel)
%         Lines_ifit{end+1}=sprintf('p.Substrate%i = ''RH''',i);
        %Lines_instr.COATING{end+1}=sprintf('\tmValues%d%s = %s ;',i,panel);
        Lines_instr.DECLARE{end+1}=sprintf('double mValues%d%s;',i,panel);
        Declared{end+1}=sprintf('mValues%d%s',i,panel);
        Lines_instr.DECLARE{end+1}=sprintf('double elementLength%d;',i);
        Declared{end+1}=sprintf('elementLength%d',i);

        if useMSteps==1
            Lines_instr.COATING{end+1}=sprintf('    mValues%i%s = round(mValues%i%s*%2.4f)/%2.4f ;',i,panel,i,panel,1/coatingOptions.deltaM,1/coatingOptions.deltaM);
            Lines_instr.COATING{end+1}=sprintf('    elementLength%d = length%d/Nsegments%i ;',i,i,i);
        else
            Lines_instr.COATING{end+1}=sprintf('    elementLength%d = length%d/%d ;',i,i,segments);
        end
        
        if coatingOptions.remove_m_under > 0
            Lines_instr.COATING{end+1}=sprintf('\tif (mValues%i%s < %d)',i,panel,coatingOptions.remove_m_under);
            Lines_instr.COATING{end+1}=sprintf('\t\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\tmValues%i%s = %d;',i,panel,coatingOptions.remove_m_under);
            Lines_instr.COATING{end+1}=sprintf('\t\t}');
        end

        if coatingOptions.remove_m_over > 0
            Lines_instr.COATING{end+1}=sprintf('\tif (mValues%i%s > %d)',i,panel,coatingOptions.remove_m_over);
            Lines_instr.COATING{end+1}=sprintf('\t\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\tmValues%i%s = %d;',i,panel,coatingOptions.remove_m_over);
            Lines_instr.COATING{end+1}=sprintf('\t\t}');
        end
        % stepScan overwrite
            Lines_instr.COATING{end+1}=sprintf('\tif (stepScan == 1)');
            Lines_instr.COATING{end+1}=sprintf('\t\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\tif (correctionMValue%i%s[0] > 0.1)',i,panel);
            Lines_instr.COATING{end+1}=sprintf('\t\t{');
            Lines_instr.COATING{end+1}=sprintf('\t\t\t\tmValues%i%s = correctionMValue%i%s[0];',i,panel,i,panel);
            Lines_instr.COATING{end+1}=sprintf('\t\t}');
            Lines_instr.DECLARE{end+1}=sprintf('double correctionMValue%i%s;',i,panel);
            Lines_instr.COATING{end+1}=sprintf('\t\t}');
        if coatingOptions.makeFree==1
        else
            [Lines_instr,Lines_ifit] = coatingPrice(i,n,panel,coatingOptions,Lines_ifit,Lines_instr,numberInThisRun,freePanels)
        end
        
    end
    
    %% Substrate price
    
    if n==1
        Lines_instr.DECLARE{end+1}=sprintf('int nsegments%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double nSegments%i;',i);
        Lines_instr.DECLARE{end+1}=sprintf('double Nsegments%i;',i);
        Lines_instr.COATINGfirst{end+1}=sprintf('nsegments%i=round(length%i/%d);',i,i,coatingOptions.mStepLength);
        Lines_instr.COATINGfirst{end+1}=sprintf('nSegments%i=round(length%i/%d);',i,i,coatingOptions.mStepLength);
        Lines_instr.COATINGfirst{end+1}=sprintf('Nsegments%i=round(length%i/%d);',i,i,coatingOptions.mStepLength);
        
        
        if strcmp(coatingOptions.substrateMode,'area') || strcmp(coatingOptions.coatingMode,'area')  %%%
            [Lines_instr] = Area(i,panel,coatingOptions,Lines_instr)
        end
    end
    [Lines_instr,Lines_ifit] = SubstratePrice(i,n,panel,coatingOptions,Lines_ifit,Lines_instr)
    
    %% Summary

%     if i==1 && n==1; Lines_instr.SUMMARY{end+1}='int len = strlen(file_name);';end
%     if i==1 && n==1; Lines_instr.SUMMARY{end+1}='file_name[len-13] = ''\0'';';end
%     if i==1 && n==1; Lines_instr.SUMMARY{end+1}='file_name=strcat(file_name,"_CoatingWriter_output");';end
%     if i==1 && n==1; Lines_instr.SUMMARY{end+1}='file_name=strcat(file_name,".txt");';end
%     if i==1 && n==1; Lines_instr.SUMMARY{end+1}='FILE *f = fopen(file_name, "w");';end
    if i==1 && n==1
     Lines_instr.DECLARE{end+1}='char summaryName[150];';
     Lines_instr.SUMMARY{end+1}='sprintf(summaryName,"CoatingWriter_output%s.txt",scanname);';
     Lines_instr.SUMMARY{end+1}='FILE *f = fopen(summaryName, "w");';
     Lines_instr.SUMMARY{end+1}='fprintf(f, "Summary of the coating distributions and prices\n");';
     Lines_instr.SUMMARY{end+1}='fprintf(f, "Created with CoatingWriter (by Martin Olsen, martinolsen.mo@gmail.com)\n");';
     Lines_instr.SUMMARY{end+1}='fprintf(f, "---------------------------------\n");';
     Lines_instr.SUMMARY{end+1}='fprintf(f, "Price=%2.2f \n",sectionPrices);';
     Lines_instr.SUMMARY{end+1}='fprintf(f, "---------------------------------\n");';
    end
    if n==1 
       Lines_instr.SUMMARY{end+1}=['fprintf(f, "Segment' num2str(i) ' \n");'  ]
       Lines_instr.SUMMARY{end+1}='fprintf(f, "---------------------------------\n");';
       if coatingOptions.makeFree==1
          Lines_instr.SUMMARY{end+1}=['fprintf(f, "Note that segment ' num2str(i) ' is not incluted in price calculations\n");']; 
       end
       
    end
    if strcmp(coatingOptions.segmentType,'E') || strcmp(coatingOptions.segmentType,'P')
        if n==1
            Lines_instr.SUMMARY{end+1}=['fprintf(f, "-- Segment ' num2str(i) ' is a ' coatingOptions.segmentType '-type guide -- \n");'] 
            Lines_instr.SUMMARY{end+1}=['fprintf(f, "Section Price: %2.2f\n",segment' num2str(i) 'Price);']; 
            Lines_instr.SUMMARY{end+1}='fprintf(f, "Element length, vertical area, horizontal area, mLeft,mRight,mTop,mBottom\n");' 
            
            Lines_instr.SUMMARY{end+1}=['for (segment =1 ; segment < nSegments' num2str(i) '+1 ; segment++)' ]
            Lines_instr.SUMMARY{end+1}='{' 
            if strcmp(coatingOptions.fixSides,'HV') 
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,  %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",elementLength' num2str(i) '[segment-1],2*Segment' num2str(i) 'VerticalArea[segment-1],2*Segment' num2str(i) 'HorizontalArea[segment-1],mValues' num2str(i) 'horizontal[segment-1],mValues' num2str(i) 'horizontal[segment-1],mValues' num2str(i) 'vertical[segment-1],mValues' num2str(i) 'vertical[segment-1] );' ]
            elseif strcmp(coatingOptions.fixSides,'none')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,  %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",elementLength' num2str(i) '[segment-1],2*Segment' num2str(i) 'VerticalArea[segment-1],2*Segment' num2str(i) 'HorizontalArea[segment-1],mValues' num2str(i) 'left[segment-1],mValues' num2str(i) 'right[segment-1],mValues' num2str(i) 'top[segment-1],mValues' num2str(i) 'bottom[segment-1] );' ]
            elseif strcmp(coatingOptions.fixSides,'all')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,  %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",elementLength' num2str(i) '[segment-1],2*Segment' num2str(i) 'VerticalArea[segment-1],2*Segment' num2str(i) 'HorizontalArea[segment-1],mValues' num2str(i) '[segment-1],mValues' num2str(i) '[segment-1],mValues' num2str(i) '[segment-1],mValues' num2str(i) '[segment-1] );' ]
            elseif strcmp(coatingOptions.fixSides,'curve')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,  %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",elementLength' num2str(i) '[segment-1],2*Segment' num2str(i) 'VerticalArea[segment-1],2*Segment' num2str(i) 'HorizontalArea[segment-1],mValues' num2str(i) 'left[segment-1],mValues' num2str(i) 'right[segment-1],mValues' num2str(i) 'vertical[segment-1],mValues' num2str(i) 'vertical[segment-1] );' ]
            end
            
            Lines_instr.SUMMARY{end+1}='}' 
        end
    elseif strcmp(coatingOptions.segmentType,'S') 
        if n==1
            Lines_instr.SUMMARY{end+1}=['fprintf(f, "-- Segment ' num2str(i) ' is a ' coatingOptions.segmentType '-type guide -- \n");'] 
            Lines_instr.SUMMARY{end+1}=['fprintf(f, "Section Price: %2.2f\n",segment' num2str(i) 'Price);']; 
            Lines_instr.SUMMARY{end+1}='fprintf(f, "Element length, vertical area, horizontal area, mLeft,mRight,mTop,mBottom\n");' 
            if strcmp(coatingOptions.fixSides,'HV') 
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) 'horizontal,mValues' num2str(i) 'horizontal,mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical );' ]
            elseif strcmp(coatingOptions.fixSides,'none')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) 'left,mValues' num2str(i) 'right,mValues' num2str(i) 'top,mValues' num2str(i) 'bottom );' ]
            elseif strcmp(coatingOptions.fixSides,'all')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ' );' ]
            elseif strcmp(coatingOptions.fixSides,'curve')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) 'left,mValues' num2str(i) 'right,mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical );' ]
            end
        end
            elseif strcmp(coatingOptions.segmentType,'C') 
        if n==1
            Lines_instr.SUMMARY{end+1}=['fprintf(f, "-- Segment ' num2str(i) ' is a ' coatingOptions.segmentType '-type guide -- \n");'] 
            Lines_instr.SUMMARY{end+1}=['fprintf(f, "Section Price: %2.2f\n",segment' num2str(i) 'Price);']; 
            Lines_instr.SUMMARY{end+1}='fprintf(f, "Element length, vertical area, horizontal area,mInside,mOutside,mTop,mBottom\n");' 
            if strcmp(coatingOptions.fixSides,'HV') 
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) 'horizontal,mValues' num2str(i) 'horizontal,mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical );' ]
            elseif strcmp(coatingOptions.fixSides,'none')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) 'left,mValues' num2str(i) 'right,mValues' num2str(i) 'top,mValues' num2str(i) 'bottom );' ]
            elseif strcmp(coatingOptions.fixSides,'all')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f ,   %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ' );' ]
            elseif strcmp(coatingOptions.fixSides,'curve')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f, "%2.2f, %2.4f, %2.4f  ,  %2.2f,   %2.2f,   %2.2f,   %2.2f, \n",length' num2str(i) ',2*Segment' num2str(i) 'VerticalArea,2*Segment' num2str(i) 'HorizontalArea,mValues' num2str(i) 'inside,mValues' num2str(i) 'outside,mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical );' ]
            end
        end
    end
   if i==coatingOptions.nInput && n==freePanels
        Lines_instr.SUMMARY{end+1}='fclose(f);' 
   end
    
   
   %% Write data to "All Data File"
   
   if i==1 && n==1
     Lines_instr.SUMMARY{end+1}='if (stepScan==1);';
     Lines_instr.SUMMARY{end+1}='{';
     Lines_instr.DECLARE{end+1}='char summaryName[150];';
     Lines_instr.SUMMARY{end+1}='sprintf(summaryName,"CoatingWriter_rawData%s.txt",scanname);';
     Lines_instr.SUMMARY{end+1}='FILE *f_allData = fopen(summaryName, "w");';
     
   end
   
   
   if strcmp(coatingOptions.segmentType,'E') || strcmp(coatingOptions.segmentType,'P')
        if n==1
            Lines_instr.SUMMARY{end+1}=['for (segment =1 ; segment < nSegments' num2str(i) '+1 ; segment++)' ]
            Lines_instr.SUMMARY{end+1}='{' 
            if strcmp(coatingOptions.fixSides,'HV') 
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',elementLength' num2str(i) '[segment-1],Segment' num2str(i) '_h1[segment-1],Segment' num2str(i) '_h2[segment-1],Segment' num2str(i) '_w1[segment-1],Segment' num2str(i) '_w2[segment-1],mValues' num2str(i) 'vertical[segment-1],mValues' num2str(i) 'vertical[segment-1],mValues' num2str(i) 'horizontal[segment-1],mValues' num2str(i) 'horizontal[segment-1],Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'none')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',elementLength' num2str(i) '[segment-1],Segment' num2str(i) '_h1[segment-1],Segment' num2str(i) '_h2[segment-1],Segment' num2str(i) '_w1[segment-1],Segment' num2str(i) '_w2[segment-1],mValues' num2str(i) 'top[segment-1],mValues' num2str(i) 'bottom[segment-1],mValues' num2str(i) 'left[segment-1],mValues' num2str(i) 'right[segment-1],Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'all')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',elementLength' num2str(i) '[segment-1],Segment' num2str(i) '_h1[segment-1],Segment' num2str(i) '_h2[segment-1],Segment' num2str(i) '_w1[segment-1],Segment' num2str(i) '_w2[segment-1],mValues' num2str(i) '[segment-1],mValues' num2str(i) '[segment-1],mValues' num2str(i) '[segment-1],mValues' num2str(i) '[segment-1],Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'curve')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',elementLength' num2str(i) '[segment-1],Segment' num2str(i) '_h1[segment-1],Segment' num2str(i) '_h2[segment-1],Segment' num2str(i) '_w1[segment-1],Segment' num2str(i) '_w2[segment-1],mValues' num2str(i) 'vertical[segment-1],mValues' num2str(i) 'vertical[segment-1],mValues' num2str(i) 'left[segment-1],mValues' num2str(i) 'right[segment-1],Substrate' num2str(i) ');' ]
            end
            Lines_instr.SUMMARY{end+1}='}' 
        end
    elseif strcmp(coatingOptions.segmentType,'S') 
        if n==1
            if strcmp(coatingOptions.fixSides,'HV') 
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical,mValues' num2str(i) 'horizontal,mValues' num2str(i) 'horizontal,Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'none')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) 'top,mValues' num2str(i) 'bottom,mValues' num2str(i) 'left,mValues' num2str(i) 'right,Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'all')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'curve')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical,mValues' num2str(i) 'left,mValues' num2str(i) 'right,Substrate' num2str(i) ');']
            end
        end
            elseif strcmp(coatingOptions.segmentType,'C') 
        if n==1
            if strcmp(coatingOptions.fixSides,'HV') 
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical,mValues' num2str(i) 'horizontal,mValues' num2str(i) 'horizontal,Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'none')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) 'top,mValues' num2str(i) 'bottom,mValues' num2str(i) 'left,mValues' num2str(i) 'right,Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'all')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',mValues' num2str(i) ',Substrate' num2str(i) ');' ]
            elseif strcmp(coatingOptions.fixSides,'curve')
                Lines_instr.SUMMARY{end+1}=['    fprintf(f_allData, "%i,' coatingOptions.segmentType ',%f,%f,%f,%f,%f,%f,%f,%f,%f,%s\n",' num2str(i) ',length' num2str(i) ',starty' num2str(i) ',endy' num2str(i) ',startx' num2str(i) ',endx' num2str(i) ',mValues' num2str(i) 'vertical,mValues' num2str(i) 'vertical,mValues' num2str(i) 'outside,mValues' num2str(i) 'inside,Substrate' num2str(i) ');']
            end
        end
    end
   
   

   if i==coatingOptions.nInput && n==freePanels
        
        Lines_instr.SUMMARY{end+1}='fclose(f_allData);'
        Lines_instr.SUMMARY{end+1}='}';
   end

end




%% define and declare from arrays (last)
for j=1:length(Define)
   if sum(cell2mat(strfind(Defined,Define{j})))==0
        Defined{end+1}=Define{j};
        Declared{end+1}=Define{j};
        Lines_instr.DEFINE{end+1}=sprintf('%s=1,',Define{j});
        Lines_instr.DECLARE{end+1}=sprintf('double %s;',char(Define{j}));
   end
end
 for j=1:length(Declare)
    if sum(cell2mat(strfind(Declared,Declare{j})))==0
        Declared{end+1}=Declare{j};
        Lines_instr.DECLARE{end+1}=sprintf('double %s;',char(Declare{j}));
    end       
 end

end

