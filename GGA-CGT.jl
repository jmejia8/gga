# GGA-CGT
# GROUPING GENETIC ALGORITHM WITH CONTROLLED GENE TRANSMISSION FOR THE BIN PACKING PROBLEM
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

function Solution(bins::Bins, status)
    s = sum( b -> sum(view(status.weight, b.w))^2, bins)
    number_of_bins = length(bins)
    fitness = s / (status.bin_capacity*number_of_bins)^2
    generation = status.generation
    fully_bins = sum( b ->  b.Bin_Fullness == status.bin_capacity, bins )
    highest_avaliable = minimum([ b.Bin_Fullness for b in bins ])

    Solution(bins,fitness,number_of_bins,generation,fully_bins,highest_avaliable)
end

function update_solution(individual::Solution, status)
    bins = individual.bins
    s = sum( b -> sum(view(status.weight, b.w))^2, bins)
    individual.number_of_bins = length(bins)
    individual.fitness = s / (status.bin_capacity*individual.number_of_bins)^2
    individual.generation = status.generation
    individual.fully_bins = sum( b ->  b.Bin_Fullness == status.bin_capacity, bins )
    individual.highest_avaliable = minimum([ b.Bin_Fullness for b in bins ])
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

    start::Float64
    end_time::Float64


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


    global_best_solution=Solution(),
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
        global_best_solution,
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

        i = 1
        for file_name in readlines("instances/instances.txt")
            @info "Solving instance: $file_name"
            data = readdlm(joinpath("instances", file_name); comments=true, comment_char='/')

            status.number_items, status.bin_capacity, status.best_solution = data[1:3]
            status.weight = sort(data[4:end], rev=true)
            status.total_accumulated_weight = sum(status.weight)

            # Sort_Descending_Weights
            status.ordered_weight = sortperm(status.weight, rev=true)
            status.L2, status.n_ = LowerBound(status)
            status.permutation = status.ordered_weight[(status.n_+1):end]

            GGA_CGT(status)

            status = Status(parms[nconf,:])


        end

    end

    print("\n\n\tEnd of process");
        


end

function GGA_CGT(status)
    # procedure GGA-CGT
    status.start = time();
    # Generate_Initial_Population() returns true if an optimal solution was found
    if !Generate_Initial_Population(status)
        
        for generation = 1:status.max_gen
            # Generation() returns 1 if an optimal solution was found
            if Generation(status)
                break
            end
        end

        # is_optimal_solution is 1 if an optimal solution was printed before
        status.global_best_solution = Find_Best_Solution(status);
        if !status.is_optimal_solution
            status.end_time = time();
            status.TotalTime = (status.end_time - status.start)
            WriteOutput(status);
        end
    end
end

function Find_Best_Solution(status)
    best = 1
    for i = 2:status.P_size
        if status.population[i].fitness < status.population[best].fitness
            best = i
        end
    end

    return deepcopy(status.population[best])
end

function Generate_Initial_Population(status)
    status.population = Solution[ Solution() for i=1:status.P_size ]

    for individual = status.population
        FF_n_(individual, status);
        
        if stop_criteria_is_met(individual,status)
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

FF(item,  individual::Solution, status) = FF(item,  individual.bins, status)

function FF(item,  bins::Bins, status)
    w = status.weight[item]
    i = findfirst( b-> b.Bin_Fullness + w <= status.bin_capacity, bins)

    if i == nothing
        push!(bins, Bin([item], w))
    else
        push!(bins[i].w, item)
        bins[i].Bin_Fullness += w
    end

end

function Copy_Solution(solution2, status, delete_solution2)
    solution = deepcopy(solution2)
    solution.generation = status.generation;

    return solution

end


# To apply the reproduction technique: Controlled selection and Controlled replacement.
# Output:
#   (1) when it finds a solution for which the size matches the L2 lower bound
#   (2) if more than 0.1*P_size individuals (solutions) have duplicated-fitness
#   (0) otherwise

function Generation(status)
    status.generation += 1
    # ==========================================================================
    #                 Controlled selection for crossover
    # ==========================================================================

    # Sort_Ascending_IndividualsFitness();
    sort!(status.population, lt = (a, b) -> a.fitness < b.fitness)
    P_size = status.P_size
    B_size = status.B_size
    p_c = status.p_c
    p_m = status.p_m


    m = P_size - floor(Int, P_size*B_size)
    random_individuals = shuffle(1:m)
    
    m = floor(Int, (1-p_c)*P_size)
    best_individuals = shuffle(m:P_size);

    k = 1
    h = length(best_individuals)
    i = P_size
    j = 1
   
    children = Solution[]

    while (i >= P_size - (p_c/2*P_size))
        f1 = best_individuals[h];   h-=1
        f2 = random_individuals[k]; k+=1
      
        if f2 == f1
            f1 = best_individuals[h]; h-=1
        end
      
        child1 = Gene_Level_Crossover_FFD(status.population[f1], status.population[f2], status)
      
        if stop_criteria_is_met(child1, status)
            return true
        end

        push!(children, child1)
      
        child2 = Gene_Level_Crossover_FFD(status.population[f2], status.population[f1], status)
      
        if stop_criteria_is_met(child1, status)
            return true
        end
        
        push!(children, child2)
        i-=1; j+=2
    end


    # ==========================================================================
    #               Controlled replacement for crossover
    # ==========================================================================
    k = 1
    for j = 1:round(Int, p_c/2*P_size - 1)
        status.population[random_individuals[k]] = children[j]
        k+=1
    end



    older_sols = findall( sol -> sol.generation < status.generation, status.population)

    j = k
    k = 1
    i = P_size
    while older_sols != nothing && i > P_size - (p_c/2*P_size)
        status.population[older_sols[k]] = children[j]
        k += 1
        i-=1
        j+=1
   end
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # aqui para abajo falta
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # ==========================================================================
    #                   Controlled selection for mutation 
    # ==========================================================================
    sort!(status.population, lt = (a, b) -> a.fitness < b.fitness)

    j = 1
    for i = reverse(1:round(Int,P_size - (p_m*P_size)))
  
        if j < P_size*B_size && status.generation - status.population[i].generation < status.life_span
      
            # ==================================================================
            #           Controlled replacement for mutation
            # ==================================================================
            status.population[j] = Copy_Solution(status.population[i], status, 0);
            Adaptive_Mutation_RP(status.population[j], status.k_cs, true, status);
          
            if stop_criteria_is_met(status.population[j], status)
                return true
            end
          
            j+=1;
        else
            Adaptive_Mutation_RP(status.population[i], status.k_ncs, false, status);
            
            if stop_criteria_is_met(status.population[i], status)
                return true
            end

        end
    end
    
    return false
end


# ==============================================================================
# To recombine two parent solutions producing a child solution.
# Input:
#   The positions in the population of the two parent solutions: father_1 and father_2
#   The position in the set of children of the child solution: child
# ==============================================================================

has_repeated_items(bin, used_items) = findfirst(used_items[bin.w]) != nothing

function Gene_Level_Crossover_FFD( father_1,  father_2, status)
    sort!(father_1.bins; lt = (a, b) -> b.Bin_Fullness < a.Bin_Fullness)
    sort!(father_2.bins; lt = (a, b) -> b.Bin_Fullness < a.Bin_Fullness)

    counter = min(father_1.number_of_bins, father_2.number_of_bins)

    bins = Bin[]

    used_items = zeros(Bool, status.number_items)

    for i = 1:counter
        if father_1.bins[i].Bin_Fullness >= father_1.bins[i].Bin_Fullness && !has_repeated_items(father_1.bins[i], used_items)
            used_items[father_1.bins[i].w] .= true
            push!(bins, father_1.bins[i])
        elseif !has_repeated_items(father_2.bins[i], used_items)
            used_items[father_1.bins[i].w] .= true
            push!(bins, father_1.bins[i])
        end
    end

    for item in findall(!, used_items)
        FF(item, bins, status)
    end

    return Solution(bins, status)


end

function stop_criteria_is_met(individual, status)
    if(individual.number_of_bins == status.L2)  
        status.end_time = time();
        status.global_best_solution = Copy_Solution( individual, status, 0);
        status.TotalTime = (status.end_time - status.start)
        status.is_optimal_solution = true
        @show status.TotalTime
        @show individual.number_of_bins
        @show status.generation
        @show status.best_solution
        WriteOutput(status);
        return true
    end
    return false
end

# To produce a small modification in a solution.
# Input:
#   The position in the population of the solution to mutate: individual
#   The rate of change to calculate the number of bins to eliminate: k
#   A value that indicates if the solution was cloned: is_cloned
function Adaptive_Mutation_RP(individual, k, is_cloned, status)

  # 
  #        number_bins,
  #        i,
  #        number_free_items = 0,
  #        free_items[ATTRIBUTES] = {0},
  #        ordered_BinFullness[ATTRIBUTES] = {0};
  #  node *p;

    pow(a,b) = a^b

    ordered_BinFullness = collect(1:individual.number_of_bins)
   
    if is_cloned
        Sort_Random(ordered_BinFullness)
    end

    # Sort_Ascending_BinFullness
    sort!(individual.bins; lt = (a, b) -> a.Bin_Fullness < b.Bin_Fullness)

    i = findfirst(b -> b.Bin_Fullness == status.bin_capacity,  individual.bins)
    if i == nothing
        i = 1
    end

    _p_ = 1 / k
    ε = (2.0 - i / individual.number_of_bins) / (i^_p_)
    pε = 1 - rand() / (i^_p_)
    number_bins = ceil(Int, i*ε*pε)

    # destruir  'number_bins' 
    free_items = zeros(Bool, status.number_items)
    for i = 1:number_bins
        bin = popfirst!(individual.bins)
        free_items[bin.w] .= true
    end
 
    individual.number_of_bins = length(individual.bins);
    number_bins = individual.number_of_bins;
 

    RP(individual, findall(free_items), status)

    update_solution(individual, status)

end



function swap(bin, weight, id1, id2, free_items, new_free_items, number_free_items, bin_capacity)
    p, s = bin.w[id1], bin.w[id2]

    sw = weight[p] + weight[s]

    for i = 1:number_free_items-1
        a = free_items[i]
        if weight[a] >= sw && bin.Bin_Fullness + weight[a]  - sw <= bin_capacity
            bin.Bin_Fullness += weight[a]  - sw
            bin.w[id1] = a

            deleteat!(bin.w, id2)
            push!(new_free_items, p, s)
            return true
        end
        
        b=free_items[i+1]
        
        if weight[b] >= sw && bin.Bin_Fullness + weight[b]  - sw <= bin_capacity
            bin.Bin_Fullness += weight[b]  - sw
            bin.w[id1] = b

            deleteat!(bin.w, id2)
            push!(new_free_items, p, s)
            return true
        end

        for j = i+1:number_free_items
            b = free_items[j]
            if weight[a] + weight[b] >= sw && bin.Bin_Fullness + weight[b] + weight[a]  - sw <= bin_capacity
                bin.Bin_Fullness += weight[a] + weight[b]   - sw
                bin.w[id1] = a
                bin.w[id2] = b

                push!(new_free_items, p, s)
                return true 
            end
        end
    end

    false
    
end

# To reinsert free items into an incomplete BPP solution.
function RP( individual,  F, status)


    individual.fitness = 0;
    individual.fully_bins = 0;

    number_free_items = length(F)

    weight = status.weight

    Sort_Random(F)

    new_free_items = Int[]
    update_sol_flag = false

    for bin in individual.bins
        for i = 1:length(bin.w)
            for j = i+1:length(bin.w)
                update_sol_flag = update_sol_flag || swap(bin, weight, i, j, F, new_free_items, number_free_items, status.bin_capacity)
            end
        end
    end
   

    shuffle!(new_free_items)
    
    for item = F
        FF(item, individual, status);
    end

    
end

function WriteOutput(status)
    
end


main()