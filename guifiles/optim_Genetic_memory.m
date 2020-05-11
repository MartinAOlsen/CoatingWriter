function [thisGeneration,options,History] = optim_Genetic(thisGeneration,options,History,criterias)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

MutateRate = 0.10;
DeathChance = 0.00;

%% Figure out how good each pop performed
for j= 1:options.populationsize
    critlist(j) = eval(['History(options.generation-1).criteria{' num2str(j) '}{1};']) * options.weights(j);
    if length(criterias) > 1
        critlist(j) = critlist(j) + eval(['History(options.generation-1).criteria{' num2str(j) '}{2};']);
        crit2(j) = eval(['History(options.generation-1).criteria{' num2str(j) '}{2};']);
    end
end




[best,bestIndex]=sort(critlist);

%% Save good solutions into a list of possible parrents
if length(History) == 2
    options.Parrents.pars = thisGeneration.pars(:,:);
    options.Parrents.crit = critlist(:);
    if length(criterias) > 1
        options.Parrents.crit2 = crit2(:);
    end
else
    tmpPars = [options.Parrents.pars, thisGeneration.pars];
    tmpCrit = [options.Parrents.crit ;critlist'];
    [b,bId] = sort(tmpCrit);
    options.Parrents.pars = tmpPars;
    options.Parrents.crit = tmpCrit;
    if length(criterias) > 1
        tmpCrit2 = [options.Parrents.crit2 ;crit2'];
        options.Parrents.crit2 = tmpCrit2;
    end
end





%% Select parrents of new population:
for i = 1 : options.populationsize
    weighted_critlist = options.Parrents.crit * options.weights(i) + options.Parrents.crit2;
    [best,bestIndex]=sort(weighted_critlist);
    
    bestIndex=flipud(bestIndex);
    bestIndex = bestIndex(1:options.populationsize);
    
    parents(1) = bestIndex(randi(numel(bestIndex)))  ;
    bestIndex_short = bestIndex(bestIndex~=parents(1));
    parents(2) = bestIndex_short(randi(numel(bestIndex_short)));
    % determine what parent gives what genome:
    for j = 1 : numel(History(1).Pars{1})
        if rand < DeathChance 
            newPop(j,i) = rand;
        else
            genomeParent = randi(2);
            newPop(j,i) = options.Parrents.pars(j,parents(genomeParent));
            if rand < MutateRate % mutate
                newPop(j,i) = rand;
            end
        end
    end
end

thisGeneration.pars = newPop;







end

