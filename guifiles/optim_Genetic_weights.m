function [thisGeneration,options,History] = optim_Genetic(thisGeneration,options,History,criterias)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Figure out how good each pop performed
for j= 1:options.populationsize
    critlist(j) = eval(['History(options.generation-1).criteria{' num2str(j) '}{1};']) * options.weights(j);
    if length(criterias) > 1
        critlist(j) = critlist(j) + eval(['History(options.generation-1).criteria{' num2str(j) '}{2};']);
        crit2(j) = eval(['History(options.generation-1).criteria{' num2str(j) '}{2};']);
    end
end




[best,bestIndex]=sort(critlist);

parentsID = bestIndex(floor(length(bestIndex)/2) : end);

%% Select parrents of new population:
for i = 1 : options.populationsize
    weighted_critlist = critlist + options.weights(i) * crit2;
    [best,bestIndex]=sort(weighted_critlist);
    
    bestIndex = bestIndex(1:4);
    
    parents(1) = i;%bestIndex(randi(numel(bestIndex)))  ;
    bestIndex_short = bestIndex(bestIndex~=parents(1));
    parents(2) = bestIndex_short(randi(numel(bestIndex_short)));
    % determine what parent gives what genome:
    for j = 1 : numel(History(1).Pars{1})
        genomeParent = randi(2);
        newPop(j,i) = thisGeneration.pars(j,parents(genomeParent)) + 0.025 * randn;
        
    end
end

thisGeneration.pars = newPop;







end

