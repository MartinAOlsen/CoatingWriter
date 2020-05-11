function[Lines_out_ifit,Declare,Define,Initialize,Coating,CoatingFirst]=CoatingWriter(coatingOptions);
addpath(genpath(cd))
overalltimer=tic
ScanType=coatingOptions.scanType;

potentialPower=1.5;
try
   potentialPower=coatingOptions.potentialPower;
catch
end
clc;
Summary={}

clear nInput;
tmpNames=fieldnames(coatingOptions);

printToGuideBot=0;
for i=1:length(tmpNames)
    if strcmp(tmpNames(i),'filePath')
        printToGuideBot=1;
        fprintf('CoatingWriter will automatically write the coating distributions in:\n%s\n',coatingOptions.filePath)
    end
end


nInput=0;
for i=1:length(tmpNames)
    if strfind(tmpNames{i},'Input')
        nInput=nInput+1;
        coatingOptions.nInput=nInput;
    end
end







%% TEMPORARY FIX TO GUIDE_BOT STANDARD 

tmp=coatingOptions;
tmp=rmfield(tmp,'scanType')
for i=1:nInput
    eval(sprintf('coatingOptions.Input%i=tmp.Input%i;',i,nInput-i+1));
end

%% Price Function Values
% This has been seperated into a seperate file
% All 
coatingOptions = PriceFunctionValues(coatingOptions)



%% Load Default values
coatingDefaults=DefaultCoatingOptions;



%% Add chosen values to lists of defaults
fprintf('Segment 1 is closest to sample and segment %i is closest to moderator\n\n',coatingOptions.nInput)
for i=1:coatingOptions.nInput
    tmp=eval(sprintf('coatingOptions.Input%i;',i));


    
%     if findstr('poly',tmp.distribution)
%         polyType=tmp.distribution(5:end);
%         if str2num(polyType)>0
%             tmp.distribution='power'
%             tmp.minNpoly=polyType
%             tmp.maxNpoly=polyType
%         end
%     end
    tmpNames = fieldnames(tmp);
    %%%% If no distribution selected give a default:
    if length(strmatch('distribution',tmpNames)) == 0; tmp.distribution=DefaultDistribution(tmp.type);end
    
    tmpDefaults=eval(sprintf('coatingDefaults.default_%s_%s',tmp.type,tmp.distribution));
  
    tmpNames = fieldnames(tmp); %% Update tmpNames
    tmpDefaults_names = fieldnames(tmpDefaults);
    
    %% Make overwirte default with input
    for j=1:length(tmpNames)
        eval(sprintf('tmpDefaults.%s=tmp.%s',char(tmpNames{j}),char(tmpNames{j})))
    end
    %default mode settings:
    fNames=fieldnames(coatingDefaults.mode);
    for j=1:length(fNames)
        eval(sprintf('tmpDefaults.%s=coatingDefaults.mode.%s',fNames{j},fNames{j}))
    end
     %default options settings:
    fNames=fieldnames(coatingDefaults.options);
    fNames_options=fieldnames(coatingOptions);
    for j=1:length(fNames)
        comp=[];
        for l=1:length(fNames_options)
            comp(l)=strcmp(fNames_options(l),fNames(j));
        end
        if sum(comp)<1;
            eval(sprintf('coatingOptions.%s=coatingDefaults.options.%s',fNames{j},fNames{j}))
        end
    end

    %% fix the price-input
%     try
%         coatingOptions.price=coatingOptions.Price
%     end
%     if length(strmatch('price',tmpNames)) > 0 || length(strmatch('Price',tmpNames)) > 0
%         if isa(tmp.price,'double') && length(tmp.price)>1
%             guess_tmp=(max(tmp.price)+min(tmp.price))/2.
%             coatingOptions.price{i}=[min(tmp.price),guess_tmp,max(tmp.price)];
%             coatingOptions.totalPrice{i}=[min(tmp.price),guess_tmp,max(tmp.price)];
%         else
%             coatingOptions.price=(tmp.price);
%             coatingOptions.totalPrice=(tmp.price);
%         end
%     elseif length(strmatch('price',tmpNames)) == 0
%         if strcmp(tmp.type,'G') == 1
%             coatingOptions.price{i}='0';
%             coatingOptions.totalPrice{i}='0';
%         else
%             coatingOptions.price={};
%             coatingOptions.price{i}={'Price'};
%         end
%     end
    %% Add chosen values to lists of defaults
    for j=1:length(tmpDefaults_names)
        tmp_defname=tmpDefaults_names(j);
        if strmatch(tmp_defname,tmpNames,'exact') > 0
            eval(sprintf('tmpDefaults.%s=tmp.%s;',tmpDefaults_names{j},tmpDefaults_names{j}))
        end
    end
    try tmpDefaults.substrateMode=coatingOptions.substrateMode; end
    try tmpDefaults.substrateArealDependancy=coatingOptions.substrateArealDependancy;end
    try tmpDefaults.coatingMode=coatingOptions.coatingMode;end
    try tmpDefaults.coatingArealDependancy=coatingOptions.coatingArealDependancy;end
%     try coatingOptions.budgetMaxError=coatingOptions.budgetMaxError;end
%     try coatingOptions.punishment=coatingOptions.punishment;end
%     try coatingOptions.potentialPower=coatingOptions.potentialPower;end
%     %modes:
%     try coatingOptions.OptimizeMode=coatingOptions.OptimizeMode;end
%     try coatingOptions.AnalyzeMode=coatingOptions.AnalyzeMode;end
    

%     if strcmp(tmpDefaults.fixSides,'all')
        eval(sprintf('Input.segment%i=tmpDefaults;',i));
         if strcmp(ScanType,'budget')  
            eval(sprintf('Input.segment%i.totalPrice=coatingOptions.Price;',i));
         end
%     elseif strcmp(tmpDefaults.fixSides,'horizontal')
%         eval(sprintf('Input.segment%i_Horizontal=tmpDefaults;',i));
%         eval(sprintf('Input.segment%i_Top=tmpDefaults;',i));
%         eval(sprintf('Input.segment%i_Bottom=tmpDefaults;',i));
%     elseif strcmp(tmpDefaults.fixSides,'vertical')
%         eval(sprintf('Input.segment%i_Vertical=tmpDefaults;',i));   
%         eval(sprintf('Input.segment%i_Left=tmpDefaults;',i));
%         eval(sprintf('Input.segment%i_Right=tmpDefaults;',i));
%     elseif strcmp(tmpDefaults.fixSides,'none')   
%         eval(sprintf('Input.segment%i_Left=tmpDefaults;',i));
%         eval(sprintf('Input.segment%i_Right=tmpDefaults;',i));
%         eval(sprintf('Input.segment%i_Top=tmpDefaults;',i));
%         eval(sprintf('Input.segment%i_Bottom=tmpDefaults;',i));
%     end
    clear tmp; clear tmpDefaults; clear tmpNames; clear tmpDefaults_names; clear tmp_defname;  %% cleanup
end

%%count how many shared prices there are
nFreePrices=0;
% for j=1:length(fieldnames(Input))
%     if strcmp((char(coatingOptions.price(j))),'price')  %% Removed "char" might cause problems later.
%         coatingOptions.nFreePrices=nFreePrices+1
%     end
% end





%% CoatingWriter settings
coatingOptions.print_to_console=0;
coatingOptions.make_debug_plot=0;  %% Bugged 
coatingOptions.save_txt_file=0;    %% Not implemented yet

% Type='E';
Lines_out_ifit={};
Declare={};
Define={};
Initialize={};
Coating={};
Monitor={};
Importance={};
Importance2={};







%% Define coating distributions and run correct CoatingLineWriter script

for i=1:nInput%length(fieldnames(Input))
    lines_ifit={};
    lines_instr={};
    ThisSegment=eval(sprintf('Input.segment%i',i));

    

    
    if strcmp(ScanType,'budget') || strcmp(ScanType,'pricescan') && ThisSegment.makeFree==0
%% ScanType=Budget
        
        Settings=eval(sprintf('Input.segment%i;',i));
        Settings.instrumentLength=coatingOptions.instrumentLength;
        Settings.price=coatingOptions.Price;
        Settings.nInput=coatingOptions.nInput;
        Settings.bisect=1;
        %% Guess mode:
        if sum(strfind(coatingOptions.OptimizeMode,'+C'),strfind(coatingOptions.OptimizeMode,'+c'))>0
            Settings.guessType='m6';
            Settings.guess.a=6;
            Settings.guess.center=0.5;
            Settings.guess.b=0;
            Settings.guess.c=0;
            Settings.guess.m=6;
            Settings.guess.Power=1;
            Settings.guess.mValues=6;
        else
            Settings.guessType='mean';
        end
        
        
                %% Set alpha fit (Jonas' fit hardcoded):
	%Settings.alpha=sprintf('%2.3f+%2.3f*m',2.1,0.24)
        if sum(strmatch('deltaM',fieldnames(eval(sprintf('coatingOptions.Input%i',i)))))>=1;
            Settings.deltaM=eval(sprintf('coatingOptions.Input%i.deltaM',i));
        end
        if sum(strmatch('mStepLength',fieldnames(eval(sprintf('coatingOptions.Input%i',i)))))>=1;
            Settings.mStepLength=eval(sprintf('coatingOptions.Input%i.mStepLength',i));
        end
       
        
        if strcmp(Settings.distribution,'none')
            fprintf('Making coating parameters for segment %i (%s) No coating is applied to this segment\n',i,ThisSegment.segmentType)
        elseif strcmp(Settings.distribution,'exp')
             Settings.formula=sprintf('A+C*exp(B*fabs(x-center))');
            if strcmp(Settings.mode,'limit')
                Settings.minA=sprintf('%d-C%i',Settings.minM,i);
                Settings.maxA=sprintf('%d-C%i',Settings.minM,i);
%                  Settings.maxB= sprintf('log(-(-maxM%i + A%i) / C%i) / abs(center%i - LH)',i,i,i,i);
%                  Settings.minB= sprintf('log(-(-minM%i + A%i) / C%i) / abs(center%i - LH)',i,i,i,i);
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
            elseif strcmp(Settings.mode,'specific')
                Settings.formula=sprintf('A+B*(x-center)*(x-center)');
                Settings.minA=sprintf('%d-C%i',Settings.minM,i);
                Settings.maxA=sprintf('%d-C%i',Settings.minM,i);
%                 Settings.maxB= Settings.maxB;
%                 Settings.minB= Settings.minB;
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
                
            elseif strcmp(Settings.mode,'eco')
                error('eco not compatible with exponential distribuion in this release')
            end

        elseif strcmp(Settings.distribution,'power') || strcmp(Settings.distribution,'polyn')
             %Budget Poly2
             Settings.formula=sprintf('B*((1/(0.5+fabs(center-0.5)))*pow((fabs(x-center)),power)) + A');
            if strcmp(Settings.mode,'limit')
                Settings.minA=sprintf('%d',Settings.minM,i);
                Settings.maxA=sprintf('%d',Settings.minM,i);
               if  Settings.minCenter== Settings.maxCenter
                    Settings.maxB= sprintf('(-A%d + %d) / abs(%d - 1) ^ 2',i, Settings.maxM, Settings.center);
                    Settings.minB= Settings.maxB;
               else
                    Settings.maxB= sprintf('(-A%d + %d) / abs(center%d - 1) ^ 2',i, Settings.maxM,i);
                    Settings.minB= Settings.maxB;
               end
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
            elseif strcmp(Settings.mode,'specific') 
%                 fNames=eval(sprintf('fieldnames(coatingOptions.Input%i)',i));
%                 for fName=1:length(fNames)
%                     eval(sprintf('Settings.%s=coatingOptions.Input%i.%s',fNames(i),i,(fNames)))
%                 end
                %Settings.remove_m_under=1;
                %Settings.remove_m_over=7;
                [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
            elseif strcmp(Settings.mode,'eco')
                error('eco not compatible with exponential distribuion in this release')
            end
            
            
            
        elseif strcmp(Settings.distribution,'linear')

%% Constant
        elseif strcmp(Settings.distribution,'constant')
            
            [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);




        end
    

    elseif strcmp(ScanType,'value') || strcmp(ScanType,'coating')   || ThisSegment.makeFree==1 || strcmp(ScanType,'manual') 
    %% ScanType=Value
    
    
        Settings=eval(sprintf('Input.segment%i;',i));
        Settings.nInput=coatingOptions.nInput;
         Settings.instrumentLength=coatingOptions.instrumentLength;
        Settings.bisect=1;
        if sum(strfind(coatingOptions.OptimizeMode,'+C'),strfind(coatingOptions.OptimizeMode,'+c'))>0
            Settings.guessType='m6';
            Settings.guess.a=6;
            Settings.guess.center=0.5;
            Settings.guess.b=0;
            Settings.guess.c=0;
            Settings.guess.m=6;
            Settings.guess.Power=1;
            Settings.guess.mValues=6;
        else
            Settings.guessType='mean';
        end
        
                %% Set alpha fit (Jonas' fit hardcoded):
	%Settings.alpha=sprintf('%2.3f+%2.3f*m',2.1,0.24)
        if sum(strmatch('deltaM',fieldnames(eval(sprintf('coatingOptions.Input%i',i)))))>=1;
            Settings.deltaM=eval(sprintf('coatingOptions.Input%i.deltaM',i));
        end
        if sum(strmatch('mStepLength',fieldnames(eval(sprintf('coatingOptions.Input%i',i)))))>=1;
            Settings.mStepLength=eval(sprintf('coatingOptions.Input%i.mStepLength',i));
        end
       
        
        if strcmp(Settings.distribution,'none')
            fprintf('Making coating parameters for segment %i (%s) No coating is applied to this segment\n',i,ThisSegment.segmentType)
        elseif strcmp(Settings.distribution,'exp')
             Settings.formula=sprintf('A+C*exp(B*fabs(x-center))');
            if strcmp(Settings.mode,'limit')
                Settings.minA=sprintf('%d-C%i',Settings.minM,i);
                Settings.maxA=sprintf('%d-C%i',Settings.minM,i);
%                  Settings.maxB= sprintf('log(-(-maxM%i + A%i) / C%i) / abs(center%i - LH)',i,i,i,i);
%                  Settings.minB= sprintf('log(-(-minM%i + A%i) / C%i) / abs(center%i - LH)',i,i,i,i);
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Value(Settings,i);
            elseif strcmp(Settings.mode,'specific')
                Settings.formula=sprintf('A+B*(x-center)*(x-center)');
                Settings.minA=sprintf('%d-C%i',Settings.minM,i);
                Settings.maxA=sprintf('%d-C%i',Settings.minM,i);
%                 Settings.maxB= Settings.maxB;
%                 Settings.minB= Settings.minB;
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Value(Settings,i);
                
            elseif strcmp(Settings.mode,'eco')
                error('eco not compatible with exponential distribuion in this release')
            end

        elseif strcmp(Settings.distribution,'power') || strcmp(Settings.distribution,'polyn')
             %Budget Poly2
             Settings.formula=sprintf('B*((1/(0.5+fabs(center-0.5)))*pow((fabs(x-center)),power)) + A');
            if strcmp(Settings.mode,'limit')
                Settings.minA=sprintf('%d',Settings.minM,i);
                Settings.maxA=sprintf('%d',Settings.minM,i);
               if  Settings.minCenter== Settings.maxCenter
                    Settings.maxB= sprintf('(-A%d + %d) / abs(%d - 1) ^ 2',i, Settings.maxM, Settings.center);
                    Settings.minB= Settings.maxB;
               else
                    Settings.maxB= sprintf('(-A%d + %d) / abs(center%d - 1) ^ 2',i, Settings.maxM,i);
                    Settings.minB= Settings.maxB;
               end
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Value(Settings,i);
            elseif strcmp(Settings.mode,'specific') 
                
                %Settings.remove_m_under=1;
                %Settings.remove_m_over=7;
                [lines_ifit,lines_instr]=CoatingLineWriter_Value(Settings,i);
            elseif strcmp(Settings.mode,'eco')
                error('eco not compatible with exponential distribuion in this release')
            end
            
            
            
        elseif strcmp(Settings.distribution,'linear')

%% Constant
        elseif strcmp(Settings.distribution,'constant')
            
            [lines_ifit,lines_instr]=CoatingLineWriter_Value(Settings,i);
            
            
% %             fprintf('Making coating parameters for segment %i (%s) with simplified CoatingWriter (for constant mvalues)\n',i,ThisSegment.segmentType)
% %             if strcmp(ThisSegment.fixSides,'all')
% %                 freePanels=1;
% %             elseif strcmp(ThisSegment.fixSides,'horizontal') || strcmp(ThisSegment.fixSides,'Horizontal') || strcmp(ThisSegment.fixSides,'horisontal')
% %                 freePanels=3;
% %                 freePanelList=[34,1,2];
% %             elseif strcmp(ThisSegment.fixSides,'vertical')
% %                 freePanels=3;
% %                 freePanelList=[12,3,4];
% %             elseif strcmp(ThisSegment.fixSides,'none')
% %                 freePanels=4;
% %                 freePanelList=[1,2,3,4];
% %             elseif strcmp(ThisSegment.fixSides,'curve') || strcmp(ThisSegment.fixSides,'curved')
% %                 freePanels=3;
% %                 freePanelList=[12,5,6];
% %             elseif strcmp(ThisSegment.fixSides,'HV') || strcmp(ThisSegment.fixSides,'VH')
% %                 freePanels=2;
% %                 freePanelList=[34,12];
% %             end
% %             lines_instr.COATING{1}=sprintf('basePrice += length%i * guidePrice;',i);
% %             
% %             for n=1:freePanels
% %                 if freePanels==1; panel='';npanels(n)=4;
% %                     elseif (freePanelList(n)==1);npanels(n)=1; panel='top';
% %                     elseif (freePanelList(n)==2);npanels(n)=1; panel='bottom';
% %                     elseif (freePanelList(n)==3);npanels(n)=1; panel='left';
% %                     elseif (freePanelList(n)==4);npanels(n)=1; panel='right';
% %                     elseif (freePanelList(n)==34);npanels(n)=2; panel='horizontal';
% %                     elseif (freePanelList(n)==12);npanels(n)=2; panel='vertical'; 
% %                     elseif (freePanelList(n)==5);npanels(n)=1; panel='inside';
% %                     elseif (freePanelList(n)==6);npanels(n)=1; panel='outside';
% %                 end
% %                 lines_ifit{n}=sprintf('p.m%i%s=[%d %2.2f %d];',i,panel,1,2,4);
% %                 lines_instr.INITIALIZE{n}='';
% %                 lines_instr.DEFINE{n}=sprintf('m%i%s=1,',i,panel);
% %                 lines_instr.DECLARE{n}=sprintf('double m%i%s;',i,panel);
% %                 if n==1;lines_instr.DECLARE{n}=sprintf('double elementlength%i;',i);end
% %                 lines_instr.DECLARE{n}=sprintf('double mValues%i%s;',i,panel);
% %                 lines_instr.COATING{end+1}=sprintf('mValues%i%s=m%i%s;',i,panel,i,panel);
% %                 lines_instr.COATING{end+1}=sprintf('sectionPrices += (pow(mValues%i%s,beta) * length%i)*%i/4;',i,panel,i,npanels(n));
% % 
% %                 lines_instr.Monitor{n}='';
% %             end
% %             lines_instr.DEFINE{end+1}=sprintf('segment%iSubstrateMaterialPrice=14,',i);
% %             lines_ifit{end+1}=sprintf('p.segment%iSubstrateMaterialPrice=''14''; %%Price of the substrate material pr. m^2 in kEuro. Base is 14 for boron and 25 for alumilium',i);
% %             lines_instr.COATING{end+1}=sprintf('TotalSubstratePrice+=segment%iSubstrateMaterialPrice*(2*(length%i*startx%i)+2*(length%i*starty%i));',i,i,i,i,i);
            
%             Settings.formula=sprintf('A');
%             [lines_ifit,lines_instr]=CoatingLineWriter_Value(Settings,i);
        end
    
    
    
    
    
   

    elseif strcmp(ScanType,'pricescan')  
    %% ScanType=PriceScan

    end
    if isempty(lines_instr)
        lines_instr.COATING{1}=sprintf('sectionPrices+=sumOfMValues;');
        lines_instr.COATING{end+1}=sprintf('priceToSave+=sumOfMValues;');
        lines_instr.COATING{end+1}=sprintf('segment%iPrice=sumOfMValues;',i);
        lines_instr.COATING{end+1}=sprintf('sumOfMValues=0;');
        lines_instr.DECLARE{1}=sprintf('double segment%iPrice;',i);
    else
        lines_instr.COATING{end+1}=sprintf('sectionPrices+=sumOfMValues;');
        lines_instr.COATING{end+1}=sprintf('priceToSave+=sumOfMValues;');
        lines_instr.COATING{end+1}=sprintf('segment%iPrice=sumOfMValues;',i);
        lines_instr.COATING{end+1}=sprintf('sumOfMValues=0;');
        lines_instr.DECLARE{end+1}=sprintf('double segment%iPrice;',i);
    end
%     % Try is a bad fix. Avoid doing this in gaps and kinks:
%     try
%         %% Set number of needed sections per element.
%         numSecPrElement=500;
%         if isfield(coatingOptions,'instrumentLength')
%             numSecPrElement=eval(['round(coatingOptions.instrumentLength/coatingOptions.Input' num2str(i) '.mStepLength)+1']);
%         end
%         %% Declare segments
%         lines_instr.DECLARE{end+1}=sprintf('double verticalSegmentSize [%i];',numSecPrElement);
%         lines_instr.DECLARE{end+1}=sprintf('double horizontalSegmentSize [%i];',numSecPrElement);    
%     end

    %% PRINT
    if length(lines_ifit)>0
    for P=1:length(lines_ifit)
        if isempty(Lines_out_ifit)
            Lines_out_ifit{1}=lines_ifit{P};
        else
            Lines_out_ifit{end+1}=lines_ifit{P};
        end
    end

    %% INSTR DEFINE
    for Q=1:length(lines_instr.DEFINE)
        if isempty(lines_instr.DEFINE)
            Define{1}='//DEFINE';
        else
            Define{end+1}=lines_instr.DEFINE{Q};
        end
    end
    %% INSTR DECLARE
    for Q=1:length(lines_instr.DECLARE)
        if isempty(lines_instr.DECLARE)
            Declare{1}='//DECLARE';
        else
            Declare{end+1}=lines_instr.DECLARE{Q};
        end
    end
    %% INSTR INITILIZE
    if i==1;
        Initialize{1}='//INITIALIZE';
        for Q=1:length(Importance)
            if isempty(Importance)
                Initialize{1}='//INITIALIZE';
            else
                Initialize{end+1}=Importance{Q};
            end
        end
    end
    
%     for Q=1:length(Importance2)
%         if isempty(Importance2)
%             Initialize{1}='//INITIALIZE';
%         else
%             Initialize{end+1}=Importance2{Q};
%         end
%     end
    
    for Q=1:length(lines_instr.INITIALIZE)
        if isempty(lines_instr.INITIALIZE)
            Initialize{1}='//INITIALIZE';
        else
            Initialize{end+1}=lines_instr.INITIALIZE{Q};
        end
    end
    
    
   for Q=1:length(lines_instr.SUMMARY)
        if isempty(lines_instr.SUMMARY)
            Summary{1}='//SUMMARY';
        else
            Summary{end+1}=lines_instr.SUMMARY{Q};
        end
    end

%     for Q=1:length(lines_instr.ALLDATA)
%         if isempty(lines_instr.ALLDATA)
%             Summary{1}='//ALLDATA';
%         else
%             Summary{end+1}=lines_instr.ALLDATA{Q};
%         end
%     end
    
    %% INSTR COATING
    Coating{1}='//COATING';
    CoatingFirst{1}='//COATINGFIRST';
    try 
        for Q=1:length(lines_instr.COATINGfirst)
        if isempty(lines_instr.COATINGfirst)
        else
            CoatingFirst{end+1}=lines_instr.COATINGfirst{Q};
        end
        end
    catch
    end
    
    for Q=1:length(lines_instr.COATING)
        if isempty(lines_instr.COATING);
        else
            Coating{end+1}=lines_instr.COATING{Q};
        end
    end
    Coating{end+1}=sprintf('MPI_MASTER(printf("\\n-----PRICE=%%2.3f k€----\\n",sectionPrices ));');


    
    
    
  
    
    %% INSTR Collect M Values
    if strcmp(coatingOptions.scanType,'manual')
       lines_instr.Monitor{end+1} =  'COMPONENT Div2d_sample_B = CoatingWriter_monitor(filename = "Monitor.txt", xwidth = sizeX,yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y,Price=sectionPrices+TotalSubstratePrice,min_lambda=WaveMin,max_lambda=WaveMax,background_min_lambda=0.1,background_max_lambda=1.5)';
       lines_instr.Monitor{end+1} = 'AT (0, 0,sample_dist) RELATIVE PREVIOUS';
       lines_instr.Monitor{end+1} = '';
       lines_instr.Monitor{end+1} = '';
        
    end
    
    
    for Q=1:length(lines_instr.Monitor)
        if isempty(lines_instr.COATING);
            Monitor{1}='//Monitor';
        else
            Monitor{end+1}=lines_instr.Monitor{Q};
        end
    end
    end

clear Settings;
end




%% Budget error
if strcmp(ScanType,'budget')  
    Coating{end+1}=sprintf('punishment=1;')
    if strcmp(coatingOptions.punishment,'potential')
%         Coating{end+1}=sprintf('if ((sectionPrices-totalPrice)> %i)',coatingOptions.budgetMaxError)
%         Coating{end+1}=sprintf('{')
%         Coating{end+1}=sprintf('sectionPrices+=(pow(fabs(sectionPrices-totalPrice-%2.2f)),%2.2d);',coatingOptions.budgetMaxError,potentialPower)
%         Coating{end+1}=sprintf('}')
%         Coating{end+1}=sprintf('else if ((sectionPrices-totalPrice)< -%i)',coatingOptions.budgetMaxError)
%         Coating{end+1}=sprintf('{')
%         Coating{end+1}=sprintf('sectionPrices+=(pow(fabs(sectionPrices-totalPrice+%2.2f)),%2.2d);',coatingOptions.budgetMaxError,potentialPower)
%         Coating{end+1}=sprintf('}')
% Moved to InsertCoating
% Coating{end+1}=sprintf('if (fabs(sectionPrices-totalPrice)> %2.2f)',coatingOptions.budgetMaxError)
%         Coating{end+1}=sprintf('{')
%         Coating{end+1}=sprintf('sectionPrices+=pow((fabs(sectionPrices-totalPrice)-%2.2f),%2.2f);',coatingOptions.budgetMaxError,potentialPower)
%         Coating{end+1}=sprintf('}')

    else
        Coating{end+1}=sprintf('if (fabs(sectionPrices-totalPrice)> %i)',coatingOptions.budgetMaxError);
        Coating{end+1}=sprintf('{');
        Coating{end+1}=sprintf('punishment=100000;');
        Coating{end+1}=sprintf('}');
    end
    Declare{end+1}=sprintf('double punishment;')
end

%% Remove scanname if none is needed
Initialize{end+1}=['if (strcmp(scanname,"noScanName") == 0)'];
Initialize{end+1}=['{'];
Initialize{end+1}=['    scanname="";'];
Initialize{end+1}=['}'];


 %% Analyse_Print creator
% lines_analyse=CWAnalyse(coatingOptions);
 


%% Add mode to result-structure
Lines_out_ifit{end+1}=['Result_this.mode=''' coatingOptions.scanType ''';']


%% Ending functon


if coatingOptions.print_to_console==1
    
    fprintf('write to ifit_coating:\n');
    for i=1:length(Lines_out_ifit)
        fprintf('%s\n',Lines_out_ifit{i});
    end
    fprintf('\n write to .instr file:\n');
    for i=1:length(Define)
        fprintf('%s\n',Define{i});
    end
    for i=1:length(Declare)
        fprintf('%s\n',Declare{i});
    end
    for i=1:length(Initialize)
        fprintf('%s\n',Initialize{i});
    end
    for i=1:length(Coating)
        fprintf('%s\n',Coating{i});
    end

    for i=1:length(Monitor)
        fprintf('%s\n',Monitor{i});
    end
else
    fprintf('\nConsole print is turned off, values can be found in lines_ifit and lines_instr\n');
    lines_instr;
    lines_ifit;
    
end

if printToGuideBot==1
        InsertCoating(Coating,CoatingFirst,Declare,Define,Initialize,Monitor,Lines_out_ifit,coatingOptions,Summary);
end



fprintf('\nCoatingWriter finished after %2.2f seconds with %i warnings\n',toc(overalltimer),0)
end




