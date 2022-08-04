%settings
global selection_no
selection_no = 1;
global crossover_no;
crossover_no=1;
global mutation_no;
mutation_no=2;
%initialising data structures, population and generation size
n = 500;
pop = zeros(n,31);
pop_new = zeros(n,30);
fitness = zeros(n,1);
trails = zeros(200,2,n);
crossover_prob = 0.8;
mutation_prob = 0.3;
Ngen = 30;
fitness_data = zeros(1,Ngen);

tic %to start recording time taken

%filling population
for z=1:n
    pop(z,1:30) = generate_ant;
end
%for loop for generations
for q=1:Ngen
    %adding fitness on to the end of encoding
    for v=1:n
        [pop(v,31),] = simulate_ant(dlmread('muir_world.txt',' '), pop(v,1:30));
    end
    
    pop = sortrows(pop,31); 
    fitness_data(q) = pop(end, 31); %saving score of the fittest for plotting later
    
    %elitism, save top 30%
    pop_new(1:(0.3*n),:) = pop(n-(0.3*n-1):n,1:30);
    pop_new_num = (0.3*n);
    %fill rest of new population
    while(pop_new_num<n)
        if(selection_no==1)
            choice1 = tournament_selection(pop(:,31));
            choice2 = tournament_selection(pop(:,31));
        end
        
        if(selection_no==2)
            weights = pop(:,31)/sum(pop(:,31));
            choice1 = roulette_wheel_selection(weights);
            choice2 = roulette_wheel_selection(weights);
        end
        parent1 = pop(choice1,1:30);
        parent2 = pop(choice2,1:30);
        if (rand<crossover_prob)
            [child1, child2] = k_point_crossover(parent1,parent2);
            
            if (rand<mutation_prob)
                if(mutation_no==1)
                    child1 = mutation(child1);
                end
                if(mutation_no==2)
                    child1 = flip_mutation(child1);
                end
            end

            if (rand<mutation_prob)
                if(mutation_no==1)
                    child2 = mutation(child2);
                end
                if(mutation_no==2)
                    child2 = flip_mutation(child2);
                end
            end
            pop_new_num = pop_new_num+1;
            pop_new(pop_new_num,:) = child1;
            pop_new_num = pop_new_num+1;
            pop_new(pop_new_num,:) = child2;
        else
            pop_new_num = pop_new_num+1;
            pop_new(pop_new_num,:) = parent1;
            pop_new_num = pop_new_num+1;
            pop_new(pop_new_num,:) = parent2;
        end
    end
    pop(:,1:30) = pop_new;
end

%find best fitness and trail for it
pop = sortrows(pop,31);
[best_fitness,trail] = simulate_ant(dlmread('muir_world.txt',' '), pop(n,1:30));

%plotting fitness over generations
hf = figure(1); set(hf,'Color',[1 1 1]);
hp = plot(1:Ngen,100*fitness_data/89,'r');
set(hp,'LineWidth',2);
axis([0 Ngen 0 100]); grid on;
xlabel('Generation number');
ylabel('Ant fitness [%]');
title('Ant fitness as a function of generation');

% read the John Moir Trail (world)
filename_world = 'muir_world.txt';
world_grid = dlmread(filename_world,' ');

% display the John Moir Trail (world)
world_grid = rot90(rot90(rot90(world_grid)));
xmax = size(world_grid,2);
ymax = size(world_grid,1);
hf = figure(2); set(hf,'Color',[1 1 1]);

for y=1:ymax
    for x=1:xmax
         if(world_grid(x,y) == 1)
             h1 = plot(x,y,'sk');
             hold on
         end
    end
end
grid on
% display the fittest Individual trail
for k=1:size(trail,1)
 h2 = plot(trail(k,2),33-trail(k,1),'*m');
 hold on
end
axis([1 32 1 32])
title_str = sprintf('John Muri Trail - Hero Ant fitness %d%% in %d generation ', uint8(100*best_fitness/89), Ngen);
title(title_str)
lh = legend([h1 h2],'Food cell','Ant movement');
set(lh,'Location','SouthEast');


toc %to end time recording

%function to generate ant encoding

function string_controller = generate_ant()
    string_controller = zeros(1,30);
    for l=1:3:30
        string_controller(l) = randi(4);
        string_controller(l+1) = randi([0,9],1);
        string_controller(l+2) = randi([0,9],1);
    end
end

%tournament selection

function choice = tournament_selection(fitness)
    [fitness, index] = sortrows(fitness);
    chr1 = randi([1 length(fitness)]);
    chr2 = randi([1 length(fitness)]);
    if(fitness(chr1)>fitness(chr2))
        choice = index(chr1);
    else
        choice = index(chr2);
    end
end

%roulette wheel selection

function choice = roulette_wheel_selection(weights)
  accumulation = cumsum(weights);
  p = rand();
  chosen_index = -1;
  for index = 1 : length(accumulation)
    if (accumulation(index) > p)
      chosen_index = index;
      break;
    end
  end
  choice = chosen_index;
end

%1 point crossover

function [child1, child2] = k_point_crossover(parent1,parent2)
    point = randi(30);
    child1 = horzcat(parent1(1:point-1),parent2(point:30));
    child2 = horzcat(parent2(1:point-1),parent1(point:30));
end

%flip mutation

function flippedchild = flip_mutation(child)
    point1 = randi([1, length(child)]);
    point2 = randi([1, length(child)]);
    while(point2 == point1)
        point2 = randi([1, length(child)]); %so that point 1 != point2
    end
    if(point2 < point1) %swap values so that point1 < point2
        [point1, point2] = swap(point1, point2);
    end
    tobeflipped = child(point1:point2);
    child(point1:point2) = fliplr(tobeflipped);
    flippedchild = child;
end

%swap function used by flip mutation

function [b, a] = swap(a, b) end

%mutation of +/-1 to one of the alleles

function child = mutation(child)
    index = randi(30);
    if (mod(index-1,3)==0)
        if(child(index)==4)
            child(index) = child(index) - 1;
        elseif(child(index)==1)
            child(index) = child(index) + 1;
        else
            child(index) = child(index) + randi([-1,1]);
        end
    else
        if(child(index)==9)
            child(index) = child(index) - 1;
        elseif(child(index)==0)
            child(index) = child(index) + 1;
        else
            child(index) = child(index) + randi([-1,1]);
        end
    end
end
            