function[Lines_out_ifit,Declare,Define,Initialize,Coating,CoatingFirst]=CoatingWriter(coatingOptions);
overalltimer=tic
ScanType=coatingOptions.scanType;


clc;

clear nInput;
tmpNames=fieldnames(coatingOptions);

printToGuideBot=0;
for i=1:length(tmpNames)
    if strcmp(tmpNames(i),'filePath')
        printToGuideBot=1;
        fprintf('CoatingWriter will automatically write the coating distributions in:\n%s\n',coatingOptions.filePath)
    end
end

for i=1:length(tmpNames)
    tmpNames{i}=tmpNames{i}(1:end-1);
end
coatingOptions.nInput=sum(strcmp('Input',tmpNames)); 
nInput=sum(strcmp('Input',tmpNames));
clear tmpNames;


% %% TEMPORARY FIX TO GUIDE_BOT STANDARD 
% 
% tmp=coatingOptions;
% tmp=rmfield(tmp,'scanType')
% for i=1:nInput
%     eval(sprintf('coatingOptions.Input%i=tmp.Input%i',i,nInput-i+1))
% end




%% Load Default values
coatingDefaults=DefaultCoatingOptions;



%% Add chosen values to lists of defaults
fprintf('Segment 1 is closest to sample and segment %i is closest to moderator\n\n',coatingOptions.nInput)
for i=1:coatingOptions.nInput
    tmp=eval(sprintf('coatingOptions.Input%i;',i));
    tmp.substrateMode=coatingOptions.substrateMode
    tmp.substrateArealDependancy=coatingOptions.substrateArealDependancy
    tmp.coatingMode=coatingOptions.coatingMode
    tmp.coatingArealDependancy=coatingOptions.coatingArealDependancy
    try
        tmp.price=coatingOptions.Price;
    end
    
    if findstr('poly',tmp.distribution)
        polyType=tmp.distribution(5:end);
        if str2num(polyType)>0
            tmp.distribution='polyN'
            tmp.minNpoly=polyType
            tmp.maxNpoly=polyType
        end
    end
    tmpNames = fieldnames(tmp);
    %%%% If no distribution selected give a default:
    if length(strmatch('distribution',tmpNames)) == 0; tmp.distribution=DefaultDistribution(tmp.type);end
    
    tmpDefaults=eval(sprintf('coatingDefaults.default_%s_%s',tmp.type,tmp.distribution));
  
    tmpNames = fieldnames(tmp); %% Update tmpNames
    tmpDefaults_names = fieldnames(tmpDefaults);
    

    %% fix the price-input
    try
        coatingOptions.price=coatingOptions.Price
    end
    if length(strmatch('price',tmpNames)) > 0 || length(strmatch('Price',tmpNames)) > 0
        if isa(tmp.price,'double') && length(tmp.price)>1
            guess_tmp=(max(tmp.price)+min(tmp.price))/2.
            coatingOptions.price{i}=[min(tmp.price),guess_tmp,max(tmp.price)];
            coatingOptions.totalPrice{i}=[min(tmp.price),guess_tmp,max(tmp.price)];
        else
            coatingOptions.price=(tmp.price);
            coatingOptions.totalPrice=(tmp.price);
        end
    elseif length(strmatch('price',tmpNames)) == 0
        if strcmp(tmp.type,'G') == 1
            coatingOptions.price{i}='0';
            coatingOptions.totalPrice{i}='0';
        else
            coatingOptions.price={};
            coatingOptions.price{i}={'Price'};
        end
    end
    %% Add chosen values to lists of defaults
    for j=1:length(tmpDefaults_names)
        tmp_defname=tmpDefaults_names(j);
        if strmatch(tmp_defname,tmpNames,'exact') > 0
            eval(sprintf('tmpDefaults.%s=tmp.%s;',tmpDefaults_names{j},tmpDefaults_names{j}))
        end
    end
    tmpDefaults.substrateMode=coatingOptions.substrateMode;
    tmpDefaults.substrateMode=coatingOptions.substrateMode
    tmpDefaults.substrateArealDependancy=coatingOptions.substrateArealDependancy
    tmpDefaults.coatingMode=coatingOptions.coatingMode
    tmpDefaults.coatingArealDependancy=coatingOptions.coatingArealDependancy

%     if strcmp(tmpDefaults.fixSides,'all')
        eval(sprintf('Input.segment%i=tmpDefaults;',i));
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
for j=1:length(fieldnames(Input))
    if strcmp((coatingOptions.price(j)),'price')  %% Removed "char" might cause problems later.
        coatingOptions.nFreePrices=nFreePrices+1
    end
end





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

    

    
    if strcmp(ScanType,'budget')  
%% ScanType=Budget
%         for j=1:nInput
%             if strcmp(eval(sprintf('coatingOptions.Input%i.type',j)),'G') || strcmp(eval(sprintf('coatingOptions.Input%i.type',j)),'K')
%             else
%                 if i==1
%                     Importance{1}=sprintf('//Total_Segments_Importance = segment%iPriceImportance',j)
%                 elseif i==length(fieldnames(Input))
%                     Importance{1}=sprintf('%s + segment%iPriceImportance;',Importance{1},j)
%                 else
%                     Importance{1}=sprintf('%s + segment%iPriceImportance',Importance{1},j)
%                 end
%             end
%         end
%         Importance{1}=sprintf('%s;',Importance{1})
     
        if strmatch('Price',fieldnames(coatingOptions))
            Settings.Price=coatingOptions.Price;
        elseif strmatch('price',fieldnames(coatingOptions))
            Settings.Price=coatingOptions.price;
        else
            error('Only one price imput supported in this version. use coatingOptions.Price=  to set a price for the coating on the entire guide')
        end
%         elseif eval(sprintf('strmatch(fieldnames(coatingOptions.Input%i),''price'')',i)) && strmatch(fieldnames(coatingOptions),'Price')==0
            
        % If only one price is given - tell CW to make a new variable
        if length(coatingOptions.Price)==sum(length(fieldnames(Input)))
            eval(sprintf('Input.segment%i.price=%2.2f',i,coatingOptions.Price))
        elseif length(coatingOptions.Price)==1
            eval(sprintf('Input.segment%i.totalPrice=%2.2f',i,coatingOptions.Price))
        else 
            error('Price input was not understood')
        end
        
        Settings=eval(sprintf('Input.segment%i;',i));
        Settings.nInput=coatingOptions.nInput;
        Settings.bisect=1;
                %% Set alpha fit (Jonas' fit hardcoded):
        Settings.alpha=sprintf('%2.3f+%2.3f*m',2.1,0.24)

        
        if sum(strmatch('deltaM',fieldnames(eval(sprintf('coatingOptions.Input%i',i)))))>=1;
            Settings.deltaM=eval(sprintf('coatingOptions.Input%i.deltaM',i));
        end
        if sum(strmatch('mStepLength',fieldnames(eval(sprintf('coatingOptions.Input%i',i)))))>=1;
            Settings.mStepLength=eval(sprintf('coatingOptions.Input%i.mStepLength',i));
        end
       
        
        if strcmp(Settings.distribution,'exp')
            %% Budget EXP
             Settings.formula=sprintf('A+C*exp(B*abs(x-center))');
            if strcmp(Settings.mode,'limit')
                Settings.minA=sprintf('%d-C%i',Settings.minM,i);
                Settings.maxA=sprintf('%d-C%i',Settings.minM,i);
%                 Settings.maxB= sprintf('log(-(-maxM%i + A%i) / C%i) / abs(center%i - LH)',i,i,i,i);
%                 Settings.minB= sprintf('log(-(-maxM%i + A%i) / C%i) / abs(center%i - LH)',i,i,i,i);
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

        elseif strcmp(Settings.distribution,'polyN') || strcmp(Settings.distribution,'polyn')
             %Budget Poly2
             Settings.formula=sprintf('B*(pow((1/(0.5+fabs(center-0.5)))*(fabs(x-center)),nPoly)) + A');
            if strcmp(Settings.mode,'limit')
                Settings.minA=sprintf('%d',Settings.minM,i);
                Settings.maxA=sprintf('%d',Settings.minM,i);
               if  Settings.minCenter== Settings.maxCenter
                    Settings.maxB= sprintf('(-A%d + %d) / fabs(%d - 1) ^ 2',i, Settings.maxM, Settings.center);
                    Settings.minB= Settings.maxB;
               else
                    Settings.maxB= sprintf('(-A%d + %d) / fabs(center%d - 1) ^ 2',i, Settings.maxM,i);
                    Settings.minB= Settings.maxB;
               end
                Settings.LH=1;
                [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
            elseif strcmp(Settings.mode,'specific')
                
                
                Settings.remove_m_under=1;
                Settings.remove_m_over=7;
                [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
            elseif strcmp(Settings.mode,'eco')
                error('eco not compatible with exponential distribuion in this release')
            end
            
            
         elseif strcmp(Settings.distribution,'constant')
             Settings.formula=sprintf('A');
             [lines_ifit,lines_instr]=CoatingLineWriter_Budget(Settings,i);
        elseif strcmp(Settings.distribution,'linear')




        end
    

    elseif strcmp(ScanType,'value')   
    %% ScanType=Value
    
    
        Settings=eval(sprintf('Input.segment%i;',i));
        Settings.nInput=coatingOptions.nInput;
        Settings.bisect=1;
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

        elseif strcmp(Settings.distribution,'polyN') || strcmp(Settings.distribution,'polyn')
             %Budget Poly2
             Settings.formula=sprintf('B*((1/(0.5+fabs(center-0.5)))*pow((fabs(x-center)),nPoly)) + A');
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
    Coating{end+1}=sprintf('printf("\\n-----PRICE=%%2.3f k€----\\n",sectionPrices );');


    %% INSTR Collect M Values
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
    Coating{end+1}=sprintf('if (fabs(sectionPrices-totalPrice)> %i)',coatingOptions.budgetMaxError)
    Coating{end+1}=sprintf('{')
    Coating{end+1}=sprintf('sectionPrices+=100000;')
    Coating{end+1}=sprintf('}')
end



 %% Analyse_Print creator
% lines_analyse=CWAnalyse(coatingOptions);
 




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
   InsertCoating(Coating,CoatingFirst,Declare,Define,Initialize,Lines_out_ifit,coatingOptions);
end

fprintf('\nCoatingWriter finished after %2.2f seconds with %i warnings\n',toc(overalltimer),0)
end




