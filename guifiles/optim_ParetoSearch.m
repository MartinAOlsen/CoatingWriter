function [thisGeneration,options,History] = optim_ParetoSearch(thisGeneration,options,History,criterias)
plotGerneration = 1;

nSwarms = floor(options.populationsize/10);
SwarmTargets = [3/(nSwarms):3/(nSwarms):3].^1.5;
C1 = 0.1;    % 
C2 = 0.2; % Random Factor
maxSpeed = 0.05 ;
minSpeed = 0.01;
ChanceOfDirectSearch = 0.01;
minDistBetweenSwarms = 0.05;
minDistBetweenSims = 0.01;

%% Create extra swarms without initial spawns
SwarmTargets(end) = 10; %% This one just look at performance, no price.

%% Initialize
if options.generation == 2
    options.particleHomeSwarm = randi(nSwarms,options.populationsize,1);
    options.swarmList = [1:nSwarms];
    options.Velocities(:,:) = (rand(options.populationsize,length(History(1).normPars{1})) - 0.5 ) * 0.1;
    options.gensSinceMove = zeros(length(SwarmTargets),1);
end

%% Get crit
for i = 1:length(History(options.generation-1).results)
    generationCriteria(i) = History(options.generation-1).results(i).intensity / History(options.generation-1).results(i).price;
end


%% Best point for each swarm:
for s = 1:length(SwarmTargets)
    allCrit=[];
    c=0;
    for i = 1:length(History(options.generation-1).results)
        for j = 1:options.generation-1
            allCrit(j,i) = History(j).results(i).intensity^(SwarmTargets(s)) / History(j).results(i).price;
            c=c+1;
            allPars(c,:) = History(j).normPars{i};
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

for s = 2:length(SwarmTargets)-2
    CritDist(s,1) = SwarmBestPrice(s)-SwarmBestPrice(s-1);
    CritDist(s,2) = SwarmBestPrice(s+1)-SwarmBestPrice(s);
    SwarmTargets(s) = (SwarmTargets(s-1) * CritDist(s,1) + SwarmTargets(s+1) * CritDist(s,2)) / (CritDist(s,1) + CritDist(s,2));
end


%% make points search
if options.generation == 2
    options.pointVelocity(:,:) = zeros(size(SwarmBest,2),size(SwarmBest,1));
end

C = 0;
for i=1:options.populationsize
    thisVelocity = options.pointVelocity(:,options.particleHomeSwarm(i));
    thisVelocity = thisVelocity + minSpeed;
    thisVelocity = thisVelocity .* ((length(thisVelocity))/sum(abs(thisVelocity)));
    
   
    
    thisSpread =  options.gensSinceMove(options.particleHomeSwarm(i)) /40;
    if options.generation == 2
        thisSpread = 10;
    end
    thisSpread = min(thisSpread,0.03);
    thisSpread = max(thisSpread,0.005);
    thisSpeed = 0;%0.01 * 1/((length(thisVelocity))/sum(abs(thisVelocity)));
    
    %thisVelocity = thisVelocity * ;
    
    while true
        for j = 1:length(thisGeneration.pars(:,i))
            thisGeneration.pars(j,i) = SwarmBest(options.particleHomeSwarm(i),j) + normrnd(thisSpeed*thisVelocity(j),thisSpread); %thisSpeed*(thisVelocity(j) + (thisSpread * randn()));
        end
        for k = 1:length(allPars(:,1))
            tmp = allPars(k,:);
            closestElementList(k) = sum(abs(tmp(:)-thisGeneration.pars(:,i)))/length(allPars(1,:));
        end
        closestElement = min(closestElementList);
        if sum(thisGeneration.pars(:,i) > 0) > 0 && sum(thisGeneration.pars(:,i) < 1) > 0 && closestElement > minDistBetweenSims
           break; 
        else
            thisSpread = thisSpread+0.01;
        end
        if closestElement < minDistBetweenSims
            C = C + 1;
        end
    end
    
   
    
    
    
    History(options.generation-1).probagated_forward(i) = 1;
    speedlist(i) = 1;
    posOut(i,:) = thisGeneration.pars(:,i);
end
fprintf('Moved %i points that were too close\n',C)

%% Overwrite if chosen to be direct search (line between two points)
for i = 1:length(SwarmTargets)
   if ChanceOfDirectSearch * options.gensSinceMove(i) > rand
       PointList=[];
       for j = 1:length(options.particleHomeSwarm )
           if options.particleHomeSwarm(j) == i
               PointList(end+1) = j;
           end
       end
       swarmList = options.swarmList;
       swarmList(swarmList==i)=[];
       toThisPoint = SwarmBest(swarmList(randi(length(swarmList))),:);
       fromThisPoint = SwarmBest(i,:);
       diff = toThisPoint-fromThisPoint;
       for j = 1:length(PointList)
           thisGeneration.pars(:,PointList(j)) = fromThisPoint + (j*diff./length(PointList));
           posOut(PointList(j),:) = thisGeneration.pars(:,PointList(j));
       end
       options.gensSinceMove(i) = 0;
   end
end


%% migrate particles
% for i=1:options.populationsize
%    if rand < ChanceOfFlip 
%         possibleflips = [1:length(SwarmTargets)];
%         possibleflips = possibleflips(possibleflips ~= options.particleHomeSwarm(i));
%         newSwarm = possibleflips(randi(length(possibleflips)));
%         fprintf('Moved particle %i from swarm %i to swarm %i\n',i,options.particleHomeSwarm(i),newSwarm)
%         options.particleHomeSwarm(i) = newSwarm;
%         
%    end
% end

%% Save last movement
options.PosHist(options.generation-1,:,:) = SwarmBest;


for j = 1:size(SwarmBest,1)
    ch = 0;
    init = options.PosHist(options.generation-1,j,:);
    init = init(:);
    for i = options.generation-2:-1:1
        test = options.PosHist(i,j,:);
        test = test(:);
        if test == init
        else
            options.pointVelocity(:,j) = init - test;
            options.gensSinceMove(j) = (options.generation-1) - i;
            break
        end
    end
end


    
    

%% Plots
if plotGerneration == 1
options.numPlots = 3;
subplot(1,3,1)
colorList = jet(length(SwarmTargets));

hold off
scatter(posOut(:,1),posOut(:,2),5,colorList(options.particleHomeSwarm(:),:))
hold on
% for i = 1:length(posOut(:,1))
%     quiver(posOut(i,1),posOut(i,2),options.Velocities(i,1),options.Velocities(i,2),'color',colorList(options.particleHomeSwarm(i),:))
% end
scatter(SwarmBest(:,1),SwarmBest(:,2),35,colorList(:,:),'filled')

xlabel('Par 1')
ylabel('Par 2')
axis([0,1,0,1])
pause(0.05)
subplot(1,3,2)
% for i = 1:length(posOut(:,1))
    for j = 1:length(posOut(1,:))
        scatter(ones(length(posOut(:,j)),1)*j,posOut(:,j),20,colorList(options.particleHomeSwarm(:),:),'filled')
        hold on
    end
% end

ylabel('val')
xlabel('Par')
xticks([1:length(options.fNames.variables)])
xticklabels(options.fNames.variables)
xtickangle(45)

pause(0.05)
    
end
end