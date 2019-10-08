function [thisGeneration,options,History] = optim_Walker(thisGeneration,options,History,criterias)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for i = 1:options.populationsize
            %% Determine how succesfull last iter was
            if options.generation > 2
                %% If solution is worse, move back before random walk
                
                
                
                moveForward = 0;

%Move if weighted criteria is better
                
                for k = 1 : options.generation - 1
                    critList(k) = 0;
                    for l = 1 : length(criterias)
                        if l == 1
                            critList(k) = critList(k) + History(k).criteria{i}{l} * options.weights(i);
                        else
                            critList(k) = critList(k) + History(k).criteria{i}{l};
                        end
                    end
                end
                [critbest critId] = max(critList);
                critId = max(critId);
                
                alpha = critList(end) / critList(options.lastMove(i));
                if critList(end) > critList(options.lastMove(i))%sum(critList(end,:) * options.weights(i)) > sum(critList(critId,:) * options.weights(i))
                    moveForward = 1;
                %elseif rand > alpha
                %    moveForward = 1;
                end
            
                        
                if moveForward == 0
                    for j = i : length(options.fNames.variables)
                        thisGeneration.pars(j,i) = History(critId).normPars{i}(j);
                    end
                    History(options.generation-1).probagated_forward(i) = 0;
                else
                    History(options.generation-1).probagated_forward(i) = 1;
                    options.lastMove(i) = options.generation - 1;
                end
            end
            
            %% Update pop
            for j = 1 : length(options.fNames.variables)
                moveMultiplier =0.025 %+ ( 0.002 * (options.generation - options.lastMove(i)));
               thisGeneration.pars(j,i) = thisGeneration.pars(j,i) + moveMultiplier * (randn);
            end
        end






end

