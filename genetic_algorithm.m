function [best, mincost] = genetic_algorithm( pop_n, crossover_p, mutation_p, data, state)
    % INPUTS:
    % pop_n: population size
    % crossover_p: crossover probability
    % mutation_p: mutation probability
    % data: binary dataset
    % state: 1 or 2 for two states of problem 1  

    App_n = size(data,1);   % applicants #
    RF_n = size(data,2);    % research feilds #
    
    % set weights of costs
    w1 = 1;    % for employment cost
    w2 = 20;   % for uncovering research feilds cost with more priority
    if state==2
        w2 = w2*1000;
    end
    
    % initialize population randomly
    pop = randi([0,1],pop_n,App_n);
    
    n=5;    % number of recent generations for convergence checking 
    bests = zeros(n, App_n);
    best = zeros(1, App_n);
    min_costs = 1000 - zeros(1,n);
    avg = 1000-zeros(1,n);
   
    gen = 0;    % generation #
    th = 3;
    if state==2
        th = th*100000;
    end
    
    while (1)
        gen = gen+1;
        cost = zeros(1, pop_n);
        for i = 1:pop_n
            RF_C = pop(i,:)*data;   % Research Feild Coverage
            
            % employment cost
            if state==1
                cost1 = sum(pop(i,:));  
            else
                cost1 = sum(RF_C)*1000;
            end            
            
            % uncovering research feilds cost
            cost2 = RF_n - length(find(RF_C));

            % total cost
            cost(i) = w1*cost1+w2*cost2;
        end
        
        [mincost, minidx] = min(cost);
        
        % save recent generation information
        bests(1:n-1,:)=bests(2:n,:);
        bests(n,:) = pop(minidx,:);
        
        min_costs(1:n-1)=min_costs(2:n);
        min_costs(n) = mincost;
        
        avg(1:n-1) = avg(2:n);
        avg(n) = mean(cost);
        
        % check for termination
        dif=avg(2:n)-avg(1:n-1);
        if abs(sum(dif))<th || isempty(find(~(dif(floor(n/2):n-1) > zeros(1,n-floor(n/2))), 1)) 
            [mincost, minidx] = min(min_costs);
            best = bests(minidx,:);
           % fprintf('number of generations: %d\n',gen);
            break;
        end
        
        % set probabilities for selection
        pr = cost/sum(cost);    % based on Cost
        pr = 1-pr;              % based on Fintness
        new_pop = zeros(pop_n, App_n);
        
        % Elitisims Reproduction 
        [~, sorted_idx] = sort(cost);
        new_pop(1:pop_n/50,:) = pop(sorted_idx(1:pop_n/50),:);
        new_pop_n = pop_n/50;
        
        while(1)
            new_pop_n = new_pop_n+1;
            if new_pop_n > pop_n
                break;
            end
            % Roulette Wheel
            cum_pr = cumsum(pr);    % cumulative sum
            parent1 = pop(find(rand<=cum_pr, 1, 'first'),:);
            parent2 = pop(find(rand<=cum_pr, 1, 'first'),:);

            % Cross Over with probability crossover_p
            if crossover_p > rand
                offspring = parent1;
                crossover_pos = round(App_n/2); 
                offspring(crossover_pos:end) = parent2(crossover_pos:end);

                % Mutation with probability mutation_p
                if mutation_p > rand
                    mutation_pos = randi([1,App_n]);
                    offspring(mutation_pos) = 1-offspring(mutation_pos);
                end
                
                new_pop(new_pop_n,:) = offspring;
            else
                new_pop(new_pop_n,:) = parent1;
                new_pop_n = new_pop_n+1;
                if new_pop_n > pop_n
                    break;
                end
                new_pop(new_pop_n,:) = parent2;
            end
        end
        pop = new_pop;
    end
    toc
  end
