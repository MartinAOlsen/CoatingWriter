function [thisGeneration,options,History] = optim_ParetoPSO(thisGeneration,options,History,criterias)
plotGerneration = 1;

nSwarms = floor(options.populationsize/10);
SwarmTargets = [3/(nSwarms):3/(nSwarms):3].^1.5;
C1 = 0.1;    % 
C2 = 0.2; % Random Factor
maxSpeed = 0.2 ;
minSpeed = 0.05;
ChanceOfFlip = 0.02;
minDistBetweenSwarms = 0.05;

%% Create extra swarms without initial spawns
SwarmTargets(end+1) = 0;
SwarmTargets(end+1) = 10000;

%% Initialize
if options.generation == 2
    options.particleHomeSwarm = randi(nSwarms,options.populationsize,1);
    options.Velocities(:,:) = (rand(options.populationsize,length(History(1).normPars{1})) - 0.5 ) * 0.1;
end

%% Get crit
for i = 1:length(History(options.generation-1).results)
    generationCriteria(i) = History(options.generation-1).results(i).intensity / History(options.generation-1).results(i).price;
end


%% Best point for each swarm:
for s = 1:length(SwarmTargets)
    allCrit=[];
    for i = 1:length(History(options.generation-1).results)
        for j = 1:options.generation-1
            allCrit(j,i) = History(j).results(i).intensity^(SwarmTargets(s)) / History(j).results(i).price;
        end
    end
    [m,I] = max(allCrit(:));
    [I_row, I_col] = ind2sub(size(allCrit),I);
    if s > 1
        firstChoice = [I_row, I_col];
        stop = 0;
        while stop < 1
            if sum(sum(allCrit)) == 0
                break;
            end
            [m,I] = max(allCrit(:));
            [I_row, I_col] = ind2sub(size(allCrit),I);
            larger = 0;
            for i=1:s-1
                if minDistBetweenSwarms > mean(abs(History(I_row).normPars{I_col}'- SwarmBest(i,:)))
                    allCrit(I_row, I_col) = 0;
                else 
                    larger = larger + 1;
                end
            end
            if larger == s-1
               stop = 1; 
            end
        end
    end
    SwarmBest(s,:) = History(I_row).normPars{I_col};
    % For debug:
    SwarmBestPrice(s) = History(I_row).results(I_col).price;
    SwarmBestIntensity(s) = History(I_row).results(I_col).intensity;
end



%% Update velocities
for i=1:options.populationsize
    pos = History(end-1).normPars{i};
    best = SwarmBest(options.particleHomeSwarm(i),:);
    dir = (best-pos') + mean(best-pos')*C2*randn(1,length(best));
    vel = options.Velocities(i,:);
    speed = SwarmBestIntensity(s) / History(end-1).results(i).intensity;
    speed = min(speed,maxSpeed);
    speed = max(speed,minSpeed);
    
    norm = 1/sum(dir);
    
    options.Velocities(i,:) = (options.Velocities(i,:) + (speed * dir)*C1 ) / (C1+1);
    
    
    History(options.generation-1).probagated_forward(i) = 1;
    speedlist(i) = sqrt(sum((C1 * speed * dir + randn(1,length(dir)) * C2).^2));
    
    posOut(i,:) = pos;
    
end

%% Update pos
for i=1:options.populationsize
    thisGeneration.pars(:,i) = History(end-1).normPars{i} + options.Velocities(i,:)';
   for j = 1:length(thisGeneration.pars(:,i))
       if thisGeneration.pars(j,i) > 1
           thisGeneration.pars(j,i) = 1;
       elseif thisGeneration.pars(j,i)  < 0
           thisGeneration.pars(j,i)  = 0;
       end
   end
end

%% migrate particles
for i=1:options.populationsize
   if rand < ChanceOfFlip 
        possibleflips = [1:length(SwarmTargets)];
        possibleflips = possibleflips(possibleflips ~= options.particleHomeSwarm(i));
        newSwarm = possibleflips(randi(length(possibleflips)));
        fprintf('Moved particle %i from swarm %i to swarm %i\n',i,options.particleHomeSwarm(i),newSwarm)
        options.particleHomeSwarm(i) = newSwarm;
        
   end
end

%% Plots
if plotGerneration == 1

subplot(1,2,1)
colorList = jet(length(SwarmTargets));

hold off
scatter(posOut(:,1),posOut(:,2),5,colorList(options.particleHomeSwarm(:),:))
hold on
for i = 1:length(posOut(:,1))
    quiver(posOut(i,1),posOut(i,2),options.Velocities(i,1),options.Velocities(i,2),'color',colorList(options.particleHomeSwarm(i),:))
end
scatter(SwarmBest(:,1),SwarmBest(:,2),35,colorList(:,:),'filled')

xlabel('Par 1')
ylabel('Par 2')
axis([0,1,0,1])
pause(0.05)

    
end
end