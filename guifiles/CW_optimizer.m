function Best = CW_optimizer(p,options,criterias,filename);
clc; close all; fprintf('\nSTARTING CW OPTIMIZER\n')
options.numPlots = 1;

if isfield(options,'RamDisk')
    if options.RamDisk ==1
        if isdir('/mnt/ramdisk') == 0
            system('sudo mkdir /mnt/ramdisk')
            system('sudo mount -t tmpfs -o rw,size=2G tmpfs /mnt/ramdisk')
        end
        tmp = regexp(options.dir,'/','split');
        options.dir = ['/mnt/ramdisk/' tmp{end-2} '/' tmp{end-1} '/' tmp{end}];
        cdRDISK = ['/mnt/ramdisk/' tmp{end-2} '/' tmp{end-1}];
        thisPath = cd;
        list=dir;
        for i = 3:length(list)
           copyfile([thisPath '/' list(i).name],[cdRDISK '/' list(i).name]) 
        end
    end
end

writeEachRun = 0;

options.generation=0;
options.lastMove(:) = ones(options.populationsize,1) ;

if options.plot == 1 ;
    fig=figure('Position',[0 0 1800 700]); 
    figAx = fig;
end
for i = 1 : length(criterias)
    Best.criteria{i} = 1;
end


%options.weights = exp(options.weightrange *  [-options.populationsize/2:options.populationsize/2] ./options.populationsize);
%options.weights = 0.5*power([2/options.populationsize:2/options.populationsize:2],options.weightrange);
options.weights = ([2/options.populationsize:2/options.populationsize:2]);
options.weights(1:floor(length(options.weights)/2)) = options.weights(1:floor(length(options.weights)/2)).^(1/options.weightrange);
options.weights(floor(length(options.weights)/2):length(options.weights)) = options.weights(floor(length(options.weights)/2):length(options.weights)).^(options.weightrange)

options.weights = options.weights;
%scatter([1:length(options.weights)],options.weights)


% options.weights = 1+randn(options.populationsize,1);
%options.weights = ones(options.populationsize,1)
%weightList = 2*([1:options.maxIter]/options.maxIter).^0.5;
weightList = ([1:options.maxIter]/options.maxIter);
options.weightsmodifier = 2*(rand(options.populationsize,1)-0.5);

%% Test if parallel computing toolbox is active
if license('test','Distrib_Computing_Toolbox')
    fprintf('\nParallel Computing Toolbox is installed\nMPI is turned off and parfor enabled\n')
    options.mpi = 1;
end


%% Compile instrument
    fprintf('compiling...')
    compile = tic;
    McString = ['mcrun'];
    McString = [McString ' ' filename];
    McString = [McString ' -n 1'];
    McString = [McString ' --dir ' options.dir];
%     if options.mpi > 1
        McString = [McString ' --mpi ' num2str(options.mpi) ''];
%     end
    McString = [McString ' -c'];
    %% Parameters
    fnames=fieldnames(p);
    for i=1:length(fnames)
        val = eval(['p.' fnames{i}]);
        if ischar(val)
            insertval = val;
        else
            insertval = num2str(val(2));
        end
        McString = [McString ' ' fnames{i} '=' insertval];
    end
    [supress,output] = unix(McString);
    %unix(McString)
    fprintf('   Done in %2.2f seconds\n',toc(compile))

%% Seperate constants and variables in p
fNames.all = fieldnames(p);


for i = 1:length(fNames.all)
    if ischar(eval(['p.' fNames.all{i}]))
        eval(['constants.' fNames.all{i} ' = p.' fNames.all{i} ';']);
    else
        if min(eval(['p.' fNames.all{i}])) == max(eval(['p.' fNames.all{i}]))
            eval(['constants.' fNames.all{i} ' = ''' num2str(max(eval(['p.' fNames.all{i}]))) ''';']);  %% If min and max is the same, it is not a viarable
            eval(['p.' fNames.all{i} ' = ''' num2str(max(eval(['p.' fNames.all{i}]))) ''';']);
        else
            eval(['variables.' fNames.all{i} ' = p.' fNames.all{i} ';']);
        end
    end
end
fNames.constants=fieldnames(constants);
fNames.variables=fieldnames(variables);

%% Optimizer loop
while options.generation < options.maxIter
    generationTimer = tic;
    
    options.generation = options.generation +1;
    %options.weights = [0.25 + 1.5/options.populationsize:1.5/options.populationsize:1.75]  * weightList(options.generation);
    %options.weights = 1 + weightList(options.generation) * options.weightsmodifier;
    %options.weights(options.weights<0) = 0;
    
    if options.generation == 1
        History(options.generation).probagated_forward = ones(options.populationsize,1);
    else
        History(options.generation).probagated_forward = zeros(options.populationsize,1);
    end
    
    
    %% Make first generation
    if options.generation == 1
        thisGeneration.pars = rand(length(fNames.variables),options.populationsize);
        
        
        %% Make first organism in population the guess from input
        for i = 1 : length(fNames.variables) 
            varSize = max(eval(['p.' fNames.variables{i}])) - min(eval(['p.' fNames.variables{i} ';']));
            %thisGeneration.pars(:,:) = ( eval(['variables.' fNames.variables{i} '(2)']) - min(eval(['p.' fNames.variables{i} ';']))) / varSize;
            thisGeneration.pars(i,1) = ( eval(['variables.' fNames.variables{i} '(2)']) - min(eval(['p.' fNames.variables{i} ';']))) / varSize;
            thisGeneration.pars(i,2) = 0;
            thisGeneration.pars(i,3) = 1;
           
        end
    end
    
    
    
    %% Probagate population  
    options.fNames = fNames;
    if options.generation > 1  
        %[thisGeneration,options,History] = optim_Walker(thisGeneration,options,History,criterias);
        %[thisGeneration,options,History] = optim_Genetic_memory(thisGeneration,options,History,criterias);
        %[thisGeneration,options,History] = optim_Genetic(thisGeneration,options,History,criterias);
        %[thisGeneration,options,History] = optim_random(thisGeneration,options,History,criterias);
        [thisGeneration,options,History] = optim_ParetoPSO(thisGeneration,options,History,criterias);
        %[thisGeneration,options,History] = optim_ParetoSearch(thisGeneration,options,History,criterias);
        
    end
    options.lastMove;
    
    %% make sure we don't go out of bounds
    for i = 1 : length(fNames.variables) 
        for j = 1 : options.populationsize
                if  thisGeneration.pars(i,j) < 0;  thisGeneration.pars(i,j) = 0; end
                if  thisGeneration.pars(i,j) > 1;  thisGeneration.pars(i,j) = 1; end
        end
    end
   
    
    %% Simulate
    %% All variables minimum value and interval to lists (before parfor due to transparency)
    fNames_variables=fNames.variables;
    variablesMatrix = zeros(length(fNames.variables),options.populationsize);
    for i = 1 : length(fNames.variables)
        varSize(i) = max(eval(['p.' fNames.variables{i}])) - min(eval(['p.' fNames.variables{i} ';']));
        varMin(i) = min(eval(['p.' fNames.variables{i} ';']));
        variablesMatrix(i,:) = thisGeneration.pars(i,:) * varSize(i) + varMin(i);
    end


    thisPars = struct;
    ParStrings = {};
    for i = 1 : options.populationsize
        ParStrings{i}='';
    %     thisPars(i) = constants;
    
        for j = 1 : length(fNames.variables)
            eval(['thisPars(i).' fNames.variables{j} ' = num2str(variablesMatrix(j,i));']); 
            eval(['ParStrings{i} = [ ParStrings{i} '' ' fNames.variables{j} '=' num2str(variablesMatrix(j,i)) '''];'] ); 
        end
        for j = 1 : length(fNames.constants)
            eval(['thisPars(i).' fNames.constants{j} ' = constants.' fNames.constants{j} ';']); 
            eval(['ParStrings{i} = [ ParStrings{i} '' ' fNames.constants{j} '=' eval(['constants.' fNames.constants{j}]) '''];'] ); 
        end
    end
%     out=struct(100);
    out={};
    
    parfor population = 1 : options.populationsize  %% Parfor
            %% Set parameters



            %% Run McStas
            run = tic;
            tmp = runMcStas(ParStrings{population},options,fNames,filename,population);
            if writeEachRun == 1
                fprintf('------- generation %3i (pop #%3i)  | ',options.generation,population)
                fprintf(' Done in %2.2f sec  |',toc(run))
                fprintf(' price = %6.0f k€ , intensity = %2.5f  | ',tmp{2,2},tmp{1,2})
                fprintf('\n')
            else
            end
            
            
    %         
            %% Extract data from run
            out(:,:,population) = tmp;

    end
    %% Collect all data from this gerneration
    % Save history
    for j = 1:length(out(1,1,:))
        for i = 1:length(out(:,1,1))
            eval(['History(' num2str(options.generation) ').results(j).' out{i,1,j} ' = ' num2str(out{i,2,j})  ' ;']);
        end
        History(options.generation).Pars(j) = {variablesMatrix(:,j)};
        History(options.generation).normPars(j) = {thisGeneration.pars(:,j)};
        %History(options.generation).AllParameterList(j) = constants;
        for k = 1:length(fNames_variables)
            eval(['History(options.generation).AllParameterList(j). ' fNames_variables{k} ' = thisPars(j).' fNames_variables{k} ' ;']);
        end
        for k = 1:length(fNames.constants)
            eval(['History(options.generation).AllParameterList(j). ' fNames.constants{k} ' = thisPars(j).' fNames.constants{k} ' ;']);
        end
    end
    
    %% Print generation details

    
    %% Calculate criteria   criterias
    for i = 1:length(out(1,1,:))
        %criteria(i,options.generation) = eval(['History(options.generation).results(i).' criterias{1} ' / History(options.generation).results(i).' criterias{2} ';']);
        for j = 1 : length(criterias)
            critCell{j} =  eval(['History(options.generation).results(i).' criterias{j} ';']);
        end
        criteria{i,options.generation} = critCell;
        %% Normalize criteria 
        if i == 1 && options.generation == 1  
            for j = 1 : length(criterias)
                critNorm(j) =  1/critCell{j};
                critAdd(j) = 0;
                
                if strcmp('price',criterias{j})
                    critNorm(j) = -critNorm(j);
                    critAdd(j) = 2;
                end
            end
        end
        %criteria{i,options.generation} = criteria{i,options.generation} * critNorm(:);
        for j = 1 : length(criterias)
           criteria{i,options.generation}{j}  = criteria{i,options.generation}{j} * critNorm(j) + critAdd(j); 
           
        end
        
    end
    History(options.generation).criteria = criteria(:,end);
    plotGradual= 1;
    %% Plot status
    if options.plot == 1
        clist = jet(options.populationsize);
           %figure(fig)
           subplot(1,options.numPlots,options.numPlots)
        if plotGradual== 1
            
            Cmap = jet(options.maxIter);
%             for j = 1:options.generation
                for i = 1 : options.populationsize
                    c1(i) = eval(['History(end).results(i).' criterias{1}]);
                    c2(i) = eval(['History(end).results(i).' criterias{2}]);
                end
                    %errorbar(c2,c1,eval(['History(options.generation).results(i).' criterias{1} '_error'])/2,'o','Color',clist(i,:),'MarkerSize',10,'Marker','.');
			if (j<options.generation)
                scatter(c2,c1,15,Cmap(options.generation,:),'filled')
			else
				scatter(c2,c1,15,Cmap(options.generation,:),'filled','LineWidth',0.6)
            end
            if options.generation == 1
                xlabel('Price [1000 euro]')
                ylabel('intensity []')
            end
                    hold on
                
%             end
            pause(0.15)
        else
        hold on
        for i = 1 : options.populationsize
            %scatter(History(options.generation).results(i).price,History(options.generation).results(i).intensity,15,clist(i,:),'filled');
            c1 = eval(['History(options.generation).results(i).' criterias{1}]);
            c2 = eval(['History(options.generation).results(i).' criterias{2}]);
            %errorbar(c2,c1,eval(['History(options.generation).results(i).' criterias{1} '_error'])/2,'o','Color',clist(i,:),'MarkerSize',10,'Marker','.');
            scatter(c2,c1,15,clist(i,:),'filled')
            if options.generation > 2 && History(options.generation-1).probagated_forward(i) == 1
                c=0;
                for j = options.generation-2 :-1: 1
                    if History(j).probagated_forward(i) == 1 && c==0;
                        propID = j;
                        j = 1;
                        c=1;
                    end
                end
                
                
                
                p1 = [eval(['History(propID).results(i).' criterias{2}])  ,  eval(['History(propID).results(i).' criterias{1}])];
                p2 = [eval(['History(options.generation-1).results(i).' criterias{2}])  ,  eval(['History(options.generation-1).results(i).' criterias{1}])];
                %plot(p1,p2,'LineStyle','--','color',clist(i,:));
                line([p1(1),p2(1)],[p1(2),p2(2)],'Color',clist(i,:),'LineStyle','--')
            end
            
        end
        
        if options.generation == 1
            for i = 1 : length(options.weights)
                legendString{i} = sprintf('pop %i , w = %2.2f',i,options.weights(i));
                
            end
            legend(legendString)
            set(fig,'DefaultLegendAutoUpdate','off')
        end
        
        %scatter(options.generation,out.divergence_background_intensity,'r','filled');
        xlabel(criterias{2})
        ylabel(criterias{1})
        title(['gen #' num2str(options.generation)])
        pause(0.1)
        %saveas(fig,['Frames/pic_' num2str(options.generation) '.jpg']);
%         print(['Frames/pic_' num2str(options.generation)],'-dpng');
%         pause(0.5)
        
    end
    
    %% Determine if new best
    for i = 1:size(criteria,1)
        newBest = 0;
        if options.generation == 1
            thisCriteria = criteria{i};
            if i == 1
               newBest = 1; 
            end
        else
            thisCriteria = criteria{i,end};
        end
        
        for j = 1 : 1%length(thisCriteria)
            if thisCriteria{1}*thisCriteria{2}  > Best.criteria{1}*Best.criteria{2}
                newBest = 1;
            end
        end
        
        
        if newBest == 1
            Best.criteria = thisCriteria;
            Best.parameters = constants;
            for j = 1 : length(fNames.variables)
                eval(['Best.parameters.' fNames.variables{j} ' = variablesMatrix(j,i);'])
            end
            Best.generation = options.generation;
            Best.pop = i;
            Best.results = History(end).results(i);
        end
    end
    

    end    
%% Print end-of-generation
fprintf('Moving to next generation -----------------------\n')
fprintf('GENERATION %i\n',options.generation)
fprintf('Best solution is gen %i (pop %1i)\n',Best.generation,Best.pop)
for i = 1 : length(criterias)
    if (Best.criteria{i}-1) > 0
        fprintf('Criteria (%s) change = +%2.2f %%\n',criterias{i},100*(Best.criteria{i}-1))
    else
        fprintf('Criteria (%s) change = %2.2f %%\n',criterias{i},100*(Best.criteria{i}-1))
    end
end
fprintf('Price = %2.0f k€ , intensity = %2.4f\n',Best.results.price,Best.results.intensity)
%
c=0;
for i= 1:numel(History)
    for j = 1:length(History(1).results)
        c=c+1;
        PriceList(c) = History(i).results(j).price;
        IList(c) = History(i).results(j).intensity;
    end
end
fprintf('Best intensities under price:\n')
fprintf('\tmax price: 1000 k€ , I = %2.4f\n',max(IList(PriceList<1000)))
fprintf('\tmax price: 2000 k€ , I = %2.4f\n',max(IList(PriceList<2000)))
fprintf('\tmax price: 3000 k€ , I = %2.4f\n',max(IList(PriceList<3000)))
fprintf('\tmax price: 4000 k€ , I = %2.4f\n',max(IList(PriceList<4000)))
%
fprintf('Generation took %2.0f seconds, avg %2.2f sec pr simulation',toc(generationTimer),toc(generationTimer)/options.populationsize)
fprintf('\n--------------------------------------------------\n')
%save('save_all.mat');
end


%fprintf('BEST: Intensity: %2.3f , Price: %2.0f k€ , sig2noise: %2.1f %% \n',BestOut.intensity,BestOut.price,100*BestOut.intensity/(BestOut.lambda_background_intensity))

Best = ManuallyChoosingOfOutput(History,p,options,fNames,filename);

end
