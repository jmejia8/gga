# GGA-CGT
# GROUPING GENETIC ALGORITHM WITH CONTROLLED GENE TRANSMISSION FOR THE BIN PACKING PROBLEM

import CSV
import DelimitedFiles.readdlm
import Random: shuffle, shuffle!

mutable struct Bin
  w::Array{Int} # weights
  Bin_Fullness::Float64 
end

Bin() = Bin(Float64[], 0.0)

Bins = Array{Bin} 

mutable struct Solution
    bins::Bins
    fitness::Float64
    number_of_bins::Int
    generation::Int # Saves the generation in which the solution was generated
    fully_bins::Int # Saves the number of bins in the solution that are fully at 100%
    highest_avaliable::Float64 # Saves the fullness of the bin with the highest avaliable capacity in the solution
end

function Solution(;bins=Bins[],
    fitness = Inf,
    number_of_bins = 0,
    generation = 0,
    fully_bins = 0,
    highest_avaliable = Inf)
    
    Solution(bins,fitness,number_of_bins,generation,fully_bins,highest_avaliable)
end

mutable struct Status
    is_optimal_solution::Bool
    save_bestSolution::Bool
    repeated_fitness::Bool
    max_gen::Int
    life_span::Int
    P_size::Int
    seed::Int

    conf::Int
    number_items::Int
    bin_capacity::Int
    generation::Int
    best_solution::Int
    n_::Int
    L2::Int
    bin_i::Int
    higher_weight::Int
    lighter_weight::Int
    ordered_weight::Array{Int, 1}
    weight::Array{Float64, 1}
    permutation::Array{Int, 1}
    items_auxiliary::Array{Int, 1}
    ordered_population::Array{Int, 1}
    best_individuals::Array{Int, 1}
    random_individuals::Array{Int, 1}


    p_m::Float64
    p_c::Float64
    k_ncs::Float64
    k_cs::Float64
    B_size::Float64
    TotalTime::Float64


    total_accumulated_weight::Float64
    weight1::Array{Float64, 1}
    _p_::Float64

    start::Int
    end_time::Int


    global_best_solution::Solution
    population::Array{Solution, 1}
    children::Array{Solution, 1}

    seed_emptybin::Int
    seed_permutation::Int
end

function Status(;
    is_optimal_solution=false,
    save_bestSolution=false,
    repeated_fitness=false,

    max_gen=0,
    life_span=0,
    P_size=0,
    seed=0,

    conf=0,
    number_items=0,
    bin_capacity=0,
    generation=0,
    best_solution=0,
    n_=0,
    L2=0,
    bin_i=0,
    higher_weight=0,
    lighter_weight=0,

    ordered_weight=Int[],
    weight=Float64[],
    permutation=Int[],
    items_auxiliary=Int[],
    ordered_population=Int[],
    best_individuals=Int[],
    random_individuals=Int[],


    p_m=0.0,
    p_c=0.0,
    k_ncs=0.0,
    k_cs=0.0,
    B_size=0.0,
    TotalTime=0.0,


    total_accumulated_weight=0,
    weight1=Float64[],
    _p_=0.0,

    start=0,
    end_time=0,


    status.global_best_solution=Solution(),
    population=Solution[],
    children=Solution[],

    seed_emptybin=0,
    seed_permutation=0,
)

    Status(is_optimal_solution,
        save_bestSolution,
        repeated_fitness,
        max_gen,
        life_span,
        P_size,
        seed,
        conf,
        number_items,
        bin_capacity,
        generation,
        best_solution,
        n_,
        L2,
        bin_i,
        higher_weight,
        lighter_weight,
        ordered_weight,
        weight,
        permutation,
        items_auxiliary,
        ordered_population,
        best_individuals,
        random_individuals,
        p_m,
        p_c,
        k_ncs,
        k_cs,
        B_size,
        TotalTime,
        total_accumulated_weight,
        weight1,
        _p_,
        start,
        end_time,
        status.global_best_solution,
        population,
        children,
        seed_emptybin,
        seed_permutation)
    
end

function Status(configuration)
    status = Status()

    status.conf    = configuration[1] 
    status.P_size  = configuration[2] 
    status.max_gen = configuration[3] 
    status.p_m     = configuration[4] 
    status.p_c     = configuration[5] 
    status.k_ncs   = configuration[6] 
    status.k_cs    = configuration[7] 
    status.B_size  = configuration[8] 
    status.life_span= configuration[9] 
    status.seed     = configuration[10] 
    status.save_bestSolution = configuration[11]

    seed_permutation = status.seed;
    seed_emptybin = status.seed;
    
    status.ordered_population = collect(1:status.P_size);
    status.random_individuals = collect(1:status.P_size);
    status.best_individuals = collect(1:status.P_size);

    status.is_optimal_solution = false;
    status.generation = 0;
    status.repeated_fitness = 0

    status
end

function main()

    # char  aux[300], nombreC[300], mystring[300];
    try
        mkdir("Solutions_GGA-CGT")
        mkdir("Details_GGA-CGT")
    catch
    end

    # READING EACH CONFIGURATION IN FILE "configurations.txt", CONTAINING THE PARAMETER VALUES FOR EACH EXPERIMENT
    input_Configurations = nothing
    try
        input_Configurations = readdlm("configurations.txt"; header=true)
    catch
        print("\n INVALID FILE");
        return false
    end

    parms, h = input_Configurations


    for nconf = 1:size(parms, 1)
        status = Status(parms[nconf,:]); st = status


        nameC = string("Solutions_GGA-CGT/GGA-CGT_(", st.conf, ").txt");
        open(nameC,"w+") do output
            write(output , "CONF\t|P|\tmax_gen\tn_m\tn_c\tk1(non-cloned_solutions)\tk2(cloned_solutions)\t|B|\tlife_span\tseed");
            write(output, "\n$(st.conf)\t$(st.P_size)\t$(st.max_gen)\t$(st.p_m)\tp_c\t$(st.k_ncs)\t$(st.k_cs)\t$(st.B_size)\t$(st.life_span)\t$(st.seed)");
            write(output, "\nInstancias \t L2 \t Bins \t Gen \t Time");
        end

        for file in readlines("instances/instances.txt")
            data = readdlm(joinpath("instances", file); comments=true, comment_char='/')

            status.number_items, status.bin_capacity, status.best_solution = data[1:3]
            status.weight = sort(data[3:end], rev=true)
            status.total_accumulated_weight = sum(status.weight)

            # Sort_Descending_Weights
            status.ordered_weight = sortperm(status.weight, rev=true)
            status.L2, status.n_ = LowerBound(status)
            status.permutation = status.ordered_weight[status.n_+1:end]

            Generation(status)

            return status
            # return status.population[end].bins
            # GGA_CGT(status)

        end

    end

    print("\n\n\tEnd of process");
        


end

function GGA_CGT(status)
    # procedure GGA-CGT
    status.start = time();
    # Generate_Initial_Population() returns true if an optimal solution was found
    if !Generate_Initial_Population(status)
        
        for generation = 1:max_gen
            # Generation() returns 1 if an optimal solution was found
            if Generation(status)
                break;
            end
            Find_Best_Solution(status);
        end

        # is_optimal_solution is 1 if an optimal solution was printed before
        if !status.is_optimal_solution
            status.end_time = time();
            status.TotalTime = (status.end_time - status.start)
            Find_Best_Solution(status);
            WriteOutput(status);
        end
    end
end

function Generate_Initial_Population(status)
    status.population = Solution[ Solution() for i=1:status.P_size ]

    for individual = status.population
        FF_n_(individual, status);
        
        individual.fitness /= individual.number_of_bins;
        
        if(individual.number_of_bins == status.L2)
            status.end_time = clock();
            status.global_best_solution = Copy_Solution(individual, false);
            
            status.TotalTime = status.end_time - status.start
            WriteOutput(status);
            status.is_optimal_solution = true
            return true
        end
    
    end
    
    return false
end

# To generate a random BPP solution with the ñ large items packed in separate bins.
# Input:
#   The position in the population of the new solution: individual

function FF_n_(individual, status)
    n_ = status.n_
    permutation = status.permutation

    status.bin_i = 1;
    individual.fully_bins = 0;
    individual.highest_avaliable = status.bin_capacity;

    if n_ > 0
        for i = 1:n_
            push!(individual.bins, Bin([i], status.weight[i]))

            if individual.bins[i].Bin_Fullness < individual.highest_avaliable
                individual.highest_avaliable = individual.bins[i].Bin_Fullness;
            end
        end


        Sort_Random(permutation);

        for j = permutation
            FF(j, individual, status);
        end

    else
        for j = shuffle(1:status.number_items)
            FF(j, individual, status);
        end
    end
    
    individual.number_of_bins = length(individual.bins)

    s = sum( b -> sum(view(status.weight, b.w))^2, individual.bins)
    individual.fitness = s / (status.bin_capacity*individual.number_of_bins)^2
end

# To sort the elements between the positions [k] and [n] of an array in random order
# Input:
#   The array to be randomized: random_array                                                                                                                 *
#   The initial random position: k
#    The final random position: n
function Sort_Random(random_array)
    shuffle!(random_array)
    return minimum(random_array)
end

 # To calculate the lower bound L2 of Martello and Toth and the ñ large items n_
function LowerBound(status)

    aux1 = aux2 = 0
    sjx=0.0; sj2=0.0; sj3=0.0
    jx=1; jp=0; jpp=0

    weight = status.weight
    ordered_weight = status.ordered_weight
    bin_capacity = status.bin_capacity
    number_items = status.number_items


    while weight[jx] > bin_capacity/2 && jx <= number_items
        jx+=1
    end
   
    n_ = jx
    if(jx == number_items)
        L2 = jx
        return L2, n_
    end
   
    if(jx == 0)
        if(mod(total_accumulated_weight, bin_capacity) >= 1)
            L2 = ceil(Int, total_accumulated_weight / bin_capacity)
        else
            L2 = floor(Int, total_accumulated_weight / bin_capacity)
            return L2, n_
        end
    else
        cj12 = jx
        for  i=jx:number_items
            sjx += weight[i]
        end

        jp = jx
        
        for i = 1:jx
            if weight[i] <= bin_capacity - weight[jx]
                jp = i
                break
            end
        end

        cj2 = jx - jp
        for i= jp:jx
            sj2 += weight[i]
        end

        jpp = jx
        sj3 = weight[ordered_weight[jpp]]
        ordered_weight[number_items] = number_items
        weight[number_items]=0
        
        while(weight[ordered_weight[jpp+1]]==weight[ordered_weight[jpp]])
            jpp+=1
            sj3 += weight[ordered_weight[jpp]]
        end
 
        L2 = cj12

        while true
          
            if(mod((sj3 + sj2),bin_capacity) >= 1)
                aux1 = ceil(Int, (sj3 + sj2)/bin_capacity - cj2)
            else
                aux1 = floor(Int, (sj3 + sj2)/bin_capacity - cj2)
            end

            if(L2 < (cj12 + aux1))
                L2 = cj12 + aux1
            end

            jpp+=1
            if(jpp < number_items)
                sj3 += weight[ordered_weight[jpp]]
                
                while(weight[ordered_weight[jpp+1]] == weight[ordered_weight[jpp]])
                  jpp+=1
                  sj3 += weight[ordered_weight[jpp]]
                end

                while(jp > 1 && weight[ordered_weight[jp-1]] <= bin_capacity - weight[ordered_weight[jpp]])
                    jp-=1
                    cj2+=1
                    sj2 += weight[ordered_weight[jp]]
                end
            
            end
        
            if(mod((sjx + sj2),bin_capacity) >= 1)
                aux2 = ceil(Int, (sjx + sj2) / bin_capacity - cj2 )
            else
                aux2 = floor(Int, (sjx + sj2) / bin_capacity - cj2 )
            end

            if(jpp <= number_items || (cj12 + aux2) > L2)
                break
            end
        end  
    end

    return L2, n_
end


# To insert an item into an incomplete BPP solution.
# Input:
#   An item to be inserted into the individual: item
# An incomplete chromosome where the item must be inserted: individual
#   The number of bins of the individual: total_bins
#   The first bin that could have sufficient available capacity to store the item: beginning
#   A value that indicates if it is the last item to be stored into the individual: is_last

function FF(item,  individual, status)
    w = status.weight[item]
    i = findfirst( bin-> bin.Bin_Fullness + w <= status.bin_capacity, individual.bins)

    if i == nothing
        push!(individual.bins, Bin([item], w))
    else
        push!(individual.bins[i].w, item)
        individual.bins[i].Bin_Fullness += w
    end


end


# To apply the reproduction technique: Controlled selection and Controlled replacement.
# Output:
#   (1) when it finds a solution for which the size matches the L2 lower bound
#   (2) if more than 0.1*P_size individuals (solutions) have duplicated-fitness
#   (0) otherwise

function Generation(status)
    # ==========================================================================
    #                 Controlled selection for crossover
    # ==========================================================================

    # Sort_Ascending_IndividualsFitness();
    sort!(population, lt = (a, b) -> a.fitness < b.fitness)

    m = P_size - floor(Int, P_size*B_size)
    random_individuals = shuffle(1:m)
    
    m = floor(Int, (1-p_c)*P_size)
    best_individuals = shuffle(m:P_size);

    k = 1
    h = P_size
    i = P_size
    j = 1
   
    while (i >= P_size - (p_c/2*P_size))
        f1 = best_individuals[h];   h-=1
        f2 = random_individuals[k]; k+=1
      
        if f2 == f1
            f1 = best_individuals[h]; h-=1
        end
      
        Gene_Level_Crossover_FFD(f1, f2, j, status);
        children[j].generation = generation + 1;
        children[j].fitness /= children[j].number_of_bins;
      
        if(children[j].number_of_bins == L2)  
            end_time = clock();
            status.global_best_solution = Copy_Solution( children[j], status, 0);
            status.global_best_solution.generation = generation+1;
            status.TotalTime = (end_time - start)
            WriteOutput(status);
            status.is_optimal_solution = true;
            return 1
        end
      
        Gene_Level_Crossover_FFD(f2, f1, j+1);
        children[j+1].generation = generation + 1;
        children[j+1].fitness /= children[j+1].number_of_bins;
      
        if(children[j+1].number_of_bins == L2)   
            end_time = clock();
            status.global_best_solution = Copy_Solution( children[j+1], status, 0);
            status.global_best_solution.generation = generation+1;
            status.TotalTime = (end_time - start);
            WriteOutput(status);
            status.is_optimal_solution = true;
            return 1
        end
        
        i-=1; j+=2
    end

    # ==========================================================================
    #               Controlled replacement for crossover
    # ==========================================================================
    k = 0;
    for j = 1:p_c/2*P_size - 1
        Copy_Solution(children[j], status, 1);
        k+=1
    end

    k = 0;

   for(i = P_size - 1; i > P_size - (p_c/2*P_size); i-=1, j+=1)
        while population[k].generation == generation + 1
            k+=1;
        end
        Copy_Solution(children[j], status, 1);
   end
    
    """
    /*****************************************************************************************************-
    ********************************Controlled selection for mutation**************************************
    ****************************************************************************************************-*/
    """
    Sort_Ascending_IndividualsFitness();
    # if(generation > 1 && repeated_fitness > 0.1*P_size)
    #   return (2);
    j = 0;
    for i = reverse(1:(P_size - (p_m*P_size)))
  
        if(j < P_size*B_size && generation+ 1 - population[i].generation < life_span)
      
            """
            /*****************************************************************************************************-
            **********************************Controlled replacement for mutation**********************************
            ****************************************************************************************************-*/
            """
            Copy_Solution(population[i], status, 0);
            Adaptive_Mutation_RP(j, k_cs, 1);
            population[j].generation = generation + 1;
            population[j].fitness /= population[j].number_of_bins;
          
            if(population[j].number_of_bins == L2)
                end_time = clock();
                status.global_best_solution = Copy_Solution( population[j], status, 0);
                status.global_best_solution.fitness = population[j].fitness;;
                status.global_best_solution.generation = generation+1;
                status.global_best_solution.number_of_bins = population[j].number_of_bins;
                status.global_best_solution.fully_bins = population[j].fully_bins;
                status.TotalTime = (end_time - start);
                WriteOutput(status);
                status.is_optimal_solution = true;
                return 1
            end
          
            j+=1;
        else
            Adaptive_Mutation_RP(i, k_ncs, 0);
            population[i].generation = generation + 1;
            population[i].fitness /= population[i].number_of_bins;
            if population[i].number_of_bins == L2
                end_time = clock();
                status.global_best_solution = Copy_Solution( population[i], status, 0);
                status.global_best_solution.fitness = population[i].fitness;;
                status.global_best_solution.generation = generation+1;
                status.global_best_solution.number_of_bins = population[i].number_of_bins;
                status.global_best_solution.fully_bins = population[i].fully_bins;
                status.TotalTime = (end_time - start);
                WriteOutput(status);
                status.is_optimal_solution = true;
                return 1
            end
        end
    end
    
    return 0;
end


# ==============================================================================
# To recombine two parent solutions producing a child solution.
# Input:
#   The positions in the population of the two parent solutions: father_1 and father_2
#   The position in the set of children of the child solution: child
# ==============================================================================
function Gene_Level_Crossover_FFD( father_1,  father_2,  child, status)

  #   k,
 #      counter,
 #      k2 = 0,
 #      ban = 1,
 #         items[ATTRIBUTES] = {0},
 #         free_items[ATTRIBUTES] = {0};
 #   children[child].highest_avaliable = bin_capacity;

    if(population[father_1].number_of_bins > population[father_2].number_of_bins)
        counter = population[father_1].number_of_bins;
    else
        counter = population[father_2].number_of_bins;
    end

     *random_order1 = new  [counter];
     *random_order2 = new  [counter];


    for(k = 1:counter)
        random_order1[k] = k;
        random_order2[k] = k;
    end

    Sort_Random(random_order1,0, population[father_1].number_of_bins);
    Sort_Random(random_order2,0, population[father_2].number_of_bins);
    Sort_Descending_BinFullness(random_order1, father_1);
    Sort_Descending_BinFullness(random_order2, father_2);

    for(k = 1:population[father_1].number_of_bins)
   
        if(population[father_1][random_order1[k]].Bin_Fullness >= population[father_2][random_order2[k]].Bin_Fullness)
            
            ban = Used_Items(father_1, random_order1[k], items);
            if (ban == 1)
                children[child][k2].w.clone_linked_list(population[father_1][random_order1[k]].w);
                children[child][k2+=1].Bin_Fullness = population[father_1][random_order1[k]].Bin_Fullness;
            
                if(children[child][k2-1].Bin_Fullness < children[child].highest_avaliable)
                    children[child].highest_avaliable = children[child][k2-1].Bin_Fullness;
                end
            end
            
            if(population[father_2][random_order2[k]].Bin_Fullness > 0)
                ban = Used_Items(father_2, random_order2[k], items);
                if (ban == 1)
                    children[child][k2].w.clone_linked_list(population[father_2][random_order2[k]].w);
                    children[child][k2+=1].Bin_Fullness = population[father_2][random_order2[k]].Bin_Fullness;

                    if(children[child][k2-1].Bin_Fullness < children[child].highest_avaliable)
                        children[child].highest_avaliable = children[child][k2-1].Bin_Fullness;
                    end
                end
            end
        else
            if(population[father_2][random_order2[k]].Bin_Fullness > 0)
                ban = Used_Items(father_2, random_order2[k], items);
                if (ban == 1)
                    children[child][k2].w.clone_linked_list(population[father_2][random_order2[k]].w);
                    children[child][k2+=1].Bin_Fullness = population[father_2][random_order2[k]].Bin_Fullness;
                    
                    if(children[child][k2-1].Bin_Fullness < children[child].highest_avaliable)
                        children[child].highest_avaliable = children[child][k2-1].Bin_Fullness;
                    end
                end
            end

            ban = Used_Items(father_1, random_order1[k], items);
            if (ban == 1)
                children[child][k2].w.clone_linked_list(population[father_1][random_order1[k]].w);
                children[child][k2+=1].Bin_Fullness = population[father_1][random_order1[k]].Bin_Fullness;
                if(children[child][k2-1].Bin_Fullness < children[child].highest_avaliable)
                    children[child].highest_avaliable = children[child][k2-1].Bin_Fullness;
                end
            end
        end
    end


    k = 0;
    for(counter = 0; counter < number_items; counter+=1)
        if(items[ordered_weight[counter]] == 0)
            free_items [k+=1] = ordered_weight[counter];
        end
    end
   
    if(k > 0)
        bin_i = 0;
        for(counter = 1:k-1)
            FF(free_items[counter], children[child], k2, bin_i,0);
            FF(free_items[counter], children[child], k2, bin_i,1);
        end
    else
        for(k = 0; k < k2; k+=1)
            children[child].fitness += pow((children[child][k].Bin_Fullness / bin_capacity), 2);
        end
    end
    children[child][number_items+1].Bin_Fullness = k2;
    free(random_order1);
    free(random_order2);
end


main()