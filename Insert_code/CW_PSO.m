function [parsOut,MonitorBest,EXITFLAG,OUTPUT] = CW_PSO(instrName,p,OPTIONS)
% [pars,monitor,m,o]
tic
% Seperate optimization vars and static pars:
fnameP = fieldnames(p);
co = 0;
for i = 1:length(fnameP)
    if ischar(eval(['p.' fnameP{i}]))
        eval(['OPTIONS.pars.' fnameP{i} ' = p.' fnameP{i} ';']);
    else  
        co = co +1;
        parsName(co) = fnameP(i);
        LB(co) = eval(['p.' fnameP{i} '(1);']);
        X0(co) = eval(['p.' fnameP{i} '(2);']);
        UB(co) = eval(['p.' fnameP{i} '(3);']);
    end

end


%% Set defaults options:
if isfield(OPTIONS,'PopulationSize')==0
    OPTIONS.PopulationSize = 25;
end
if isfield(OPTIONS,'MaxIter')==0
    OPTIONS.MaxIter = 5000;
end
if isfield(OPTIONS,'TolFun')==0
    OPTIONS.TolFun = 1e-8;
end
if isfield(OPTIONS,'TolX')==0
    OPTIONS.TolX = 1e-3;
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
if isfield(OPTIONS,'OPTIONS.MinFunEvals')==0
    OPTIONS.MinFunEvals = 1;
end
if isfield(OPTIONS,'OPTIONS.MinIter')==0
    OPTIONS.MinIter = 1;
end

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
        fprintf('CW_PSO starting iter %i. Number %i in generation %i',nFUN_EVALS,j,i)
        [FITNESS(j,i),mon(j)] = RunSwarm(instrName,SWARM(j,:,i),parsName,OPTIONS);
        nFUN_EVALS = nFUN_EVALS + 1;
    end
    
    % identify particle's location at which the best FITNESS has been achieved (PBEST)
    
    for j = 1:OPTIONS.PopulationSize,
        
        [PBEST_FITNESS(j,i),INDEX_PBEST(j,i)] = min(FITNESS(j,:));
        PBEST(j,:,i) = SWARM(j,:,INDEX_PBEST(j,i));
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
    
    if isempty(OPTIONS.MinFunEvals) || nFUN_EVALS >= OPTIONS.MinFunEvals
    if OPTIONS.TolX >0 & max(max(abs(diff(SWARM(:,:,i),1,1)))) < OPTIONS.TolX,
        message='Change in X less than the specified tolerance (TolX).';
        EXITFLAG = -5;
    end
    if OPTIONS.TolFun & abs(min(FITNESS(:,i))-GBEST_FITNESS(i)) < abs(OPTIONS.TolFun) ...
       & abs(min(FITNESS(:,i))-GBEST_FITNESS(i)) > 0
      EXITFLAG=-12;
      message = [ 'Termination function change tolerance criteria reached (options.TolFun=' ...
                num2str(OPTIONS.TolFun) ')' ];
    end
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
parsOut = p;
counter = 0;
for i = 1:length(fnameP)
    if ischar(eval(['p.' fnameP{i}]))
    else  
        counter= counter +1;
        eval(['parsOut.' fnameP{i} '=' num2str(X(counter)) ';']);
    end

end

fprintf('\n------------ OPTIMIZATION FINISHED ------------\n')
try
    time_str = SECS2HMS(toc);
    fprintf('Optimization ran %i iterations in %s. Avg time was %2.2f seconds pr iteration\n',nFUN_EVALS,time_str,toc/nFUN_EVALS)
end
fprintf('Optimization started at criteria %2.3f and ended at %2.3f\n',FITNESS(1,1),FVAL)
fprintf('An improvement of %2.4f%%\n',(FVAL-FITNESS(1,1))/FITNESS(1,1))
fprintf('Exiting Optimizer')
return
end

% ==============================================================================

% COST FUNCTION EVALUATION
% ------------------------



function [YTRY,mon] = RunSwarm(instrName,particle,parsName,OPTIONS);
    mcStasString = ['mcrun ' instrName ' '];
    if OPTIONS.gravitation 
        mcStasString = [mcStasString '-g '];
    end
    mcStasString = [mcStasString '-n ' num2str(OPTIONS.ncount) ' '];
    mcStasString = [mcStasString '--dir ' OPTIONS.dir ' '];
    mcStasString = [mcStasString '--mpi ' num2str(OPTIONS.mpi) ' '];
    for i = 1:length(particle)
        mcStasString = [mcStasString parsName{i} '=' num2str(particle(i)) ' '];
    end
    fnames = fieldnames(OPTIONS.pars);
    for i = 1:length(fnames)
        mcStasString = [mcStasString fnames{i} '=' eval(['OPTIONS.pars.' fnames{i}]) ' '];
    end
    
    system(mcStasString);
    
    mon = iLoad([OPTIONS.dir '/' OPTIONS.monitors '.dat']);
    
    system(['rm -rf ' OPTIONS.dir]);
    
    YTRY = -abs(sum(sum(mon.data)));

return
end

% 
% 
% function [YTRY] = CALCULATE_COST(FUN,PTRY,LB,UB,varargin)
% 
% global NDIM nFUN_EVALS
% 
% % add one to number of function evaluations
% nFUN_EVALS = nFUN_EVALS + 1;
% 
% for i = 1:NDIM,
%     % check lower bounds
%     if PTRY(i) < LB(i),
%         YTRY = 1e12+(LB(i)-PTRY(i))*1e6;
%         return
%     end
%     % check upper bounds
%     if PTRY(i) > UB(i),
%         YTRY = 1e12+(PTRY(i)-UB(i))*1e6;
%         return
%     end
% end
% 
% % calculate cost associated with PTRY
% YTRY = feval(FUN,PTRY,varargin{:});
% 
% return

