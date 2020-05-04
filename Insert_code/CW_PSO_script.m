fprintf('Starting CW_PSO Optimizer')
%% Timer
optimTimer = tic;

% Seperate optimization vars and static pars:
fnameP = fieldnames(p_optim);
co = 0;
for i = 1:length(fnameP)
    if ischar(eval(['p_optim.' fnameP{i}]))
        eval(['OPTIONS.pars.' fnameP{i} ' = p_optim.' fnameP{i} ';']);
    else  
        co = co +1;
        parsName(co) = fnameP(i);
        LB(co) = eval(['p_optim.' fnameP{i} '(1);']);
        X0(co) = eval(['p_optim.' fnameP{i} '(2);']);
        UB(co) = eval(['p_optim.' fnameP{i} '(3);']);
    end

end

if isfield(OPTIONS,'TolX') && ischar(OPTIONS.TolX) 
	OPTIONS.TolX_char = OPTIONS.TolX;
	OPTIONS = rmfield(OPTIONS,'TolX');
end
if isfield(OPTIONS,'TolFun') && ischar(OPTIONS.TolFun) 
	OPTIONS.TolFun_char = OPTIONS.TolFun;
	OPTIONS = rmfield(OPTIONS,'TolFun');
end


%% Set defaults options:
if isfield(OPTIONS,'TolX') && ischar(OPTIONS.TolX) 
	OPTIONS.TolX_char = OPTIONS.TolX;
	OPTIONS = rmfield(OPTIONS,'TolX')
end
if isfield(OPTIONS,'TolFun') && ischar(OPTIONS.TolFun) 
	OPTIONS.TolFun_char = OPTIONS.TolFun;
	OPTIONS = rmfield(OPTIONS,'TolFun')
end
if isfield(OPTIONS,'TolX_char')==0
    OPTIONS.TolX_char = '0.07%';
end
if isfield(OPTIONS,'TolFun_char')==0
    OPTIONS.TolFun_char = '0.03%';
end

if isfield(OPTIONS,'PopulationSize')==0
    %% Default pop size is 3 times dimentionality, but max is chosen to be 100 due to performance
	% Min is 25
    OPTIONS.PopulationSize = min(4*co , 100);
    OPTIONS.PopulationSize = max(OPTIONS.PopulationSize , 25);
end
if isfield(OPTIONS,'MaxIter')==0
    OPTIONS.MaxIter = 5000;
end
if isfield(OPTIONS,'MaxFunEvals')==0
    OPTIONS.MaxFunEvals = 5000;
end
if isfield(OPTIONS,'SwarmC1')==0
    OPTIONS.SwarmC1 = 2.8;
end
if isfield(OPTIONS,'SwarmC2')==0 
    OPTIONS.SwarmC2 = 1.3;
end
if isfield(OPTIONS,'MinFunEvals')==0
    OPTIONS.MinFunEvals = 150;
end
if isfield(OPTIONS,'MinIter')==0
    OPTIONS.MinIter = 10;
end
if isfield(OPTIONS,'PlotProgress')==0
    if select == 1
        OPTIONS.PlotProgress = 1;
    else
        OPTIONS.PlotProgress = 0;
    end
end

OPTIONS.MinFunEvals = max(OPTIONS.MinIter * OPTIONS.PopulationSize , OPTIONS.MinFunEvals);



%% Remove price files 
% Migh have an entry with a price from compilation.... This will remove it
% if that is the case
try; system('rm -rf priceListPunished.txt'); end
try; system('rm -rf priceList.txt'); end

nFUN_EVALS = 0;

% set EXITFLAG to default value

EXITFLAG = 0;

% determine number of variables to be optimized

NDIM = length(X0);

% seed the random number generator

rand('state',sum(100*clock));

% initialize swarm (each row of swarm corresponds to one particle)

SWARM = zeros(OPTIONS.PopulationSize,NDIM,OPTIONS.MaxIter);

for i = 1:OPTIONS.PopulationSize,
    if i == 1,
        SWARM(1,:,1) = X0(:)';
    else
        SWARM(i,:,1) = LB(:)' + rand(1,NDIM).*(UB(:)'-LB(:)');
    end
end

% initialize VELOCITIES

VELOCITIES = zeros(OPTIONS.PopulationSize,NDIM,OPTIONS.MaxIter);

% initialize FITNESS, PBEST_FITNESS, GBEST_FITNESS, INDEX_PBEST, index_gbest_particle and INDEX_GBEST_ITERATION

FITNESS = nan(OPTIONS.PopulationSize,OPTIONS.MaxIter);
PBEST = nan(OPTIONS.PopulationSize,NDIM,OPTIONS.MaxIter);
GBEST = nan(OPTIONS.MaxIter,NDIM);
PBEST_FITNESS = nan(OPTIONS.PopulationSize,OPTIONS.MaxIter);
GBEST_FITNESS = nan(OPTIONS.PopulationSize,1);
INDEX_PBEST = nan(OPTIONS.PopulationSize,OPTIONS.MaxIter);
INDEX_GBEST_PARTICLE = nan(OPTIONS.MaxIter,1);
INDEX_GBEST_ITERATION = nan(OPTIONS.MaxIter,1);

% calculate constriction factor from acceleration coefficients

if OPTIONS.SwarmC1+OPTIONS.SwarmC2 <= 4,
    % Display message
    if strcmp(OPTIONS.Dispplay,'iter') | strcmp(OPTIONS.Display,'notify')
      disp('Sum of Cognitive Acceleration Coefficient and Social Acceleration Coefficient is less then or equal to 4.')
      disp('Their values were adjusted automatically to satisfy this condition.');
      disp(' ')
    end
    % the values are adjusted so that the sum is equal to 4.1, keeping the ratio SwarmC1/SwarmC2 constant
    OPTIONS.SwarmC1 = OPTIONS.SwarmC1*4.1/(OPTIONS.SwarmC1+OPTIONS.SwarmC2);
    OPTIONS.SwarmC2 = OPTIONS.SwarmC2*4.1/(OPTIONS.SwarmC1+OPTIONS.SwarmC2);
    % calculate constriction factor
    k = 1; % k can take values between 0 and 1, but is usually set to one (Montes de Oca et al., 2006)
    OPTIONS.ConstrictionFactor = 2*k/(abs(2-(OPTIONS.SwarmC1+OPTIONS.SwarmC2)-sqrt((OPTIONS.SwarmC1+OPTIONS.SwarmC2)^2-4*(OPTIONS.SwarmC1+OPTIONS.SwarmC2))));
else
    % calculate constriction factor
    k = 1; % k can take values between 0 and 1, but is usually set to one (Montes de Oca et al., 2006)
    OPTIONS.ConstrictionFactor = 2*k/(abs(2-(OPTIONS.SwarmC1+OPTIONS.SwarmC2)-sqrt((OPTIONS.SwarmC1+OPTIONS.SwarmC2)^2-4*(OPTIONS.SwarmC1+OPTIONS.SwarmC2))));
end

% initialize counters

nITERATIONS = 0;
nFUN_EVALS = 0;
message='';

% for each iteration....

for i = 1:OPTIONS.MaxIter,
    
    % if a termination criterium was met, the value of EXITFLAG should have changed
    % from its default value of -2 to -1, 0, 1 or 2
    
    if EXITFLAG
        break
    end
    
    % calculate FITNESS values for all particles in SWARM 
    % (each row of FITNESS corresponds to the FITNESS value of one particle)
    % (each column of FITNESS corresponds to the FITNESS values of the particles in one iteration)
    clear mon
    for j = 1:OPTIONS.PopulationSize,
        fprintf('\n\n-------- CW_PSO -------- \n')
        fprintf('CW_PSO starting iter %i. Number %i in generation %i\n',nFUN_EVALS+1,j,i)
        time = toc(optimTimer);
        hours = floor(time / 3600);
        time = time - hours * 3600;
        mins = floor(time / 60);
        secs = time - mins * 60;
        fprintf('running for %1.0f hours %1.0f minutes %1.2f seconds \n',hours,mins,secs)
        fprintf('Average time for simulations: %2.2f seconds\n',toc(optimTimer)/nFUN_EVALS)
        fprintf('Improvement of %2.4f%%\n',100*(min(min(FITNESS(:,:)))-FITNESS(1,1))/FITNESS(1,1))
        fprintf('\n')
        %[FITNESS(j,i),mon(j)] = RunSwarm(instrName,SWARM(j,:,i),parsName,OPTIONS);
        %function [YTRY,mon] = RunSwarm(instrName,particle,parsName,OPTIONS);
        particle = SWARM(j,:,i);
        mcStasString = ['mcrun ' instrName ' '];
        if OPTIONS.gravitation 
            mcStasString = [mcStasString '-g '];
        end
        mcStasString = [mcStasString '-n ' num2str(OPTIONS.ncount) ' '];
        if isfield(OPTIONS,'Seed')
            mcStasString = [mcStasString '--seed ' num2str(OPTIONS.Seed) ' '];
        end
        mcStasString = [mcStasString '--dir ' OPTIONS.dir ' '];
        mcStasString = [mcStasString '--mpi ' num2str(OPTIONS.mpi) ' '];
        for I = 1:length(particle)
            mcStasString = [mcStasString parsName{I} '=' num2str(particle(I)) ' '];
        end
        fnames = fieldnames(OPTIONS.pars);
        for I = 1:length(fnames)
            mcStasString = [mcStasString fnames{I} '=' eval(['OPTIONS.pars.' fnames{I}]) ' '];
        end
        try;system(['rm -rf ' OPTIONS.dir]);end
        fprintf('\nEXCECUTING: %s\n',mcStasString)
        dos(mcStasString);

        tmp = iLoad([OPTIONS.dir '/' OPTIONS.monitors '.dat']);
        
        mon(j) = tmp;

    
        if select == 1
            FITNESS(j,i) = -abs(sum(sum(tmp.data)));
        else
            FITNESS(j,i) = -abs(sum(sum(tmp.Data.I)));
        end
	fprintf('Fitness = %2.5f\n',FITNESS(j,i))
        nFUN_EVALS = nFUN_EVALS + 1;
    end
    
    % identify particle's location at which the best FITNESS has been achieved (PBEST)
    
    for j = 1:OPTIONS.PopulationSize,
        
        [PBEST_FITNESS(j,i),INDEX_PBEST(j,i)] = min(FITNESS(j,:));
        PBEST(j,:,i) = SWARM(j,:,INDEX_PBEST(j,i));
    end
    
     % Set tol from first iteration
    if i == 1 
	if isfield(OPTIONS,'TolFun_char')
    	    OPTIONS.TolFun = abs(str2num(OPTIONS.TolFun_char(1:end-1))*FITNESS(1,1))
	end
	if isfield(OPTIONS,'TolX_char')
    	    OPTIONS.TolX = abs(str2num(OPTIONS.TolX_char(1:end-1))*(UB(:)-LB(:)))
	end
    end



    % identify the particle from the SWARM at which the best FITNESS has been achieved so far (GBEST)
        
    [GBEST_FITNESS(i),index_gbest] = min(reshape(FITNESS,numel(FITNESS),1));
    [INDEX_GBEST_PARTICLE(i),INDEX_GBEST_ITERATION(i)] = ind2sub(size(FITNESS),index_gbest);
    GBEST(i,:) = SWARM(INDEX_GBEST_PARTICLE(i),:,INDEX_GBEST_ITERATION(i));
    if INDEX_GBEST_ITERATION(i) == i
        MonitorBest = mon(INDEX_GBEST_PARTICLE(i));    
    end
    % update the VELOCITIES
    
    VELOCITIES(:,:,i+1) = OPTIONS.ConstrictionFactor.*(VELOCITIES(:,:,i) + OPTIONS.SwarmC1.*rand(OPTIONS.PopulationSize,NDIM).*(PBEST(:,:,i)-SWARM(:,:,i)) + OPTIONS.SwarmC2.*rand(OPTIONS.PopulationSize,NDIM).*(repmat(GBEST(i,:),[OPTIONS.PopulationSize 1 1 1])-SWARM(:,:,i)));
    
    % update particle positions
    
    SWARM(:,:,i+1) = SWARM(:,:,i)+VELOCITIES(:,:,i+1);
    
    % to make sure that particles stay within specified bounds...
    %   (suppose that the particle's new position is outside the boundaries,
    %    then the particle's position is adjusted by assuming that the boundary
    %    acts like a wall or mirror) (selfmade solution)
    
    for j = 1:OPTIONS.PopulationSize,
        for k = 1:NDIM,
            % check upper boundary
            if length(UB) == 1,
                if SWARM(j,k,i+1) > UB,
                    SWARM(j,k,i+1) = UB-rand*(SWARM(j,k,i+1)-UB);
                    VELOCITIES(j,k,i+1) = SWARM(j,k,i+1)-SWARM(j,k,i);
                end
            else
                if SWARM(j,k,i+1) > UB(k),
                    SWARM(j,k,i+1) = UB(k)-rand*(SWARM(j,k,i+1)-UB(k));
                    VELOCITIES(j,k,i+1) = SWARM(j,k,i+1)-SWARM(j,k,i);
                end
            end
            % check lower boundary
            if length(UB) == 1,
                if SWARM(j,k,i+1) < LB,
                    SWARM(j,k,i+1) = LB+rand*(LB-SWARM(j,k,i+1));
                    VELOCITIES(j,k,i+1) = SWARM(j,k,i+1)-SWARM(j,k,i);
                end
            else
                if SWARM(j,k,i+1) < LB(k),
                    SWARM(j,k,i+1) = LB(k)+rand*(LB(k)-SWARM(j,k,i+1));
                    VELOCITIES(j,k,i+1) = SWARM(j,k,i+1)-SWARM(j,k,i);
                end
            end
        end
    end    
    %%%%%%%%%% Print
	printStep='optimize';
	time=toc;
	try;eval(printScript);end
    %%%%%%%%%%    
    % update counters
    
    nITERATIONS = nITERATIONS+1;
    
    % give user feedback on screen if requested
    
    if 0 && strcmp(OPTIONS.Display,'iter'),
        if nITERATIONS == 1,
            disp(' Nr Iter  Nr Fun Eval    Current best function    Current worst function       Best function');
            disp(sprintf(' %5.0f     %5.0f             %12.6g              %12.6g           %15.6g',nITERATIONS,nFUN_EVALS,min(FITNESS(:,i)),max(FITNESS(:,i)),GBEST_FITNESS(i)));
        else
            disp(sprintf(' %5.0f     %5.0f             %12.6g              %12.6g           %15.6g',nITERATIONS,nFUN_EVALS,min(FITNESS(:,i)),max(FITNESS(:,i)),GBEST_FITNESS(i)));
        end
    end
    
  % end the optimization if one of the stopping criteria is met
    %% 1. difference between best and worst function evaluation in simplex is smaller than TolFun 
    %% 2. maximum difference between the coordinates of the vertices in simplex is less than TolX
    %% 3. no convergence,but maximum number of iterations has been reached
	if i > 1
        Change = abs(sum(sum((diff(GBEST_FITNESS(1:i))./GBEST_FITNESS(2:i))'.*([2:i]/i).^4)) / sum(([1:i]/i).^4));
    else
        Change = 1;
    end    


    if isempty(OPTIONS.MinFunEvals) || nFUN_EVALS >= OPTIONS.MinFunEvals
    if OPTIONS.TolX' >0 & (max(abs(diff(SWARM(:,:,i),1,1)))) < OPTIONS.TolX',
        message='Change in X less than the specified tolerance (TolX).';
        EXITFLAG = -5;
	fprintf('\nEXITFLAG: Change in X less than the specified tolerance (TolX)\n')
    end
   
    if OPTIONS.TolFun_char & Change < abs(str2num(OPTIONS.TolFun_char(1:end-1))/100)  & Change > 0
      EXITFLAG=-12;
	fprintf('\nEXITFLAG: tolerance criteria reached\n')
      message = [ 'Termination function change tolerance criteria reached (options.TolFun=' ...
                num2str(OPTIONS.TolFun) ')' ];
    end
%    if OPTIONS.TolFun & abs(min(FITNESS(:,i))-GBEST_FITNESS(i)) < abs(OPTIONS.TolFun) ...
%       & abs(min(FITNESS(:,i))-GBEST_FITNESS(i)) > 0
%      EXITFLAG=-12;
%	fprintf('\nEXITFLAG: tolerance criteria reached\n')
%      message = [ 'Termination function change tolerance criteria reached (options.TolFun=' ...
%                num2str(OPTIONS.TolFun) ')' ];
%    end
    end
    if nFUN_EVALS >= OPTIONS.MaxFunEvals
        EXITFLAG = 1;
	fprintf('\nEXITFLAG: number of evals larger than max\n')
    end
    if nITERATIONS >= OPTIONS.MaxIter
        EXITFLAG = 1;
	fprintf('\nEXITFLAG: number of iters larger than max\n')
    end
    if OPTIONS.PlotProgress == 1
	if i == 1
		close all
        	figure('Renderer', 'painters', 'Position', [10 10 900 500])
	end
        subplot(1,3,1)
        hold on
        scatter(SWARM(:,(1),i),SWARM(:,(2),i),15,'filled')

        title('pars 1 and 2')
        ylabel(parsName{1})
        xlabel(parsName{2})
        
        
        subplot(1,3,2)
        hold off
        counter = 0;
        for I = 1:i
            for J = 1:size(FITNESS,1)
                counter=counter +1;
                F(counter) = FITNESS(J,I);
            end
        end
        plot(F);
        
        m(i) = mean(FITNESS(1:OPTIONS.PopulationSize,i));
        for I = 1:i
            line([(I-1)*OPTIONS.PopulationSize,(I)*OPTIONS.PopulationSize],[m(I),m(I)],'color',[1,0,0],'linewidth',2)
        end
        
        title('Criteria History')
        xlabel('iter')
        ylabel('criteria')
        
        subplot(1,3,3)
        if i > 1
            hold on
            scatter(i,Change,25,'filled');
            line([0,i],[abs(str2num(OPTIONS.TolFun_char(1:end-1))/100),abs(str2num(OPTIONS.TolFun_char(1:end-1))/100)])
            title('Change of best')
            xlabel('change')
            ylabel('iter')
        else
            scatter(i,0,25,'filled');
        end
        pause(1)
    end

    
    if EXITFLAG
      break
    end

end

% return solution

X = GBEST(i,:);
FVAL = GBEST_FITNESS(i);

% store number of function evaluations

OUTPUT.funcCount = nFUN_EVALS;

% store number of iterations

OUTPUT.iterations = nITERATIONS;
OUTPUT.message    = message;
OUTPUT.criteriaBest = FVAL;
OUTPUT.parsBest = X; 
OUTPUT.totalTime = toc; 
OUTPUT.avgTime = toc/nFUN_EVALS; 
counter=0;
for i = 1:size(FITNESS,1)
    for j = 1:size(FITNESS,2)
        if isnan(FITNESS(i,j)) == 0
            counter=counter +1;
            OUTPUT.criteriaHistory(counter) = FITNESS(i,j);
        end
    end
end
%parsOut = p_optim;
counter = 0;
for i = 1:length(fnameP)
    if ischar(eval(['p_optim.' fnameP{i}]))
    else  
        counter= counter +1;
        eval(['parsOut.' fnameP{i} '=' num2str(X(counter)) ';']);
    end

end

%[pars,monitor,EXITFLAG,o] = [parsOut,MonitorBest,EXITFLAG,OUTPUT] 
o=OUTPUT; clear Output;
monitor=MonitorBest; clear MonitorBest;
pars=parsOut; clear parsOut;
%% Due to differences between cluster and home computer
try;monitor.Data.Parameters = monitor.Param;end
if select == 2
	for i=1:length(fnameP)
        eval(['monitor.Data.Parameters.' fnameP{i} '=p_optim.' fnameP{i} ';'])
    end
    parsFname = fieldnames(pars);
    for i=1:length(parsFname)
        eval(['monitor.Data.Parameters.' parsFname{i} '=pars.' parsFname{i} ';'])
    end
end
monitor.Data.Parameters.scanname = p.scanname;
monitor

fprintf('\n------------ OPTIMIZATION FINISHED ------------\n')
time = toc(optimTimer);
hours = floor(time / 3600);
time = time - hours * 3600;
mins = floor(time / 60);
secs = time - mins * 60;
fprintf('Optimization ran %i iterations in %1.0f hours %1.0f minutes %1.2f seconds. Avg time was %2.2f seconds pr iteration\n',nFUN_EVALS,hours,mins,secs,toc/nFUN_EVALS)
fprintf('Optimization started at criteria %2.8f and ended at %2.8f\n',FITNESS(1,1),FVAL)
fprintf('An improvement of %2.4f%%\n',100*(FVAL-FITNESS(1,1))/FITNESS(1,1))
fprintf('Exiting Optimizer')
