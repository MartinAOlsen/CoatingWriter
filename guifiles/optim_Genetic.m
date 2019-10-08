function [thisGeneration,options,History] = optim_Genetic(thisGeneration,options,History,criterias)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Figure out how good each pop performed
for j= 1:options.populationsize
    critlist(j) = eval(['History(options.generation-1).criteria{' num2str(j) '}{1};']) * options.weights(j);
    if length(criterias) > 1
        critlist(j) = critlist(j) + eval(['History(options.generation-1).criteria{' num2str(j) '}{2};']);
    end
end




[best,bestIndex]=sort(critlist);

parentsID = bestIndex(floor(length(bestIndex)/2) : end);

%% Select parrents of new population:
for i = 1 : options.populationsize
    parents(1) = parentsID(randi(numel(parentsID)))  ;
    parentsID_short = parentsID(parentsID~=parents(1));
    parents(2) = parentsID_short(randi(numel(parentsID_short)));
    % determine what parent gives what genome:
    for j = 1 : numel(History(1).Pars{1})
        genomeParent = randi(2);
        newPop(j,i) = thisGeneration.pars(j,parents(genomeParent)) + 0.01 * randn;
        
    end
    
    
    
end

thisGeneration.pars = newPop;







end

