# GGA-CGT
# GROUPING GENETIC ALGORITHM WITH CONTROLLED GENE TRANSMISSION FOR THE BIN PACKING PROBLEM
import DelimitedFiles.readdlm
import Random: shuffle, shuffle!, seed!
import Printf.@sprintf

include("io-utils.jl")
include("structures.jl")

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

        for file_name in readlines("instances/instances.txt")
            status.nameC = nameC
            @info "Solving instance: $file_name"
            status.file = file_name
            data = readdlm(joinpath("instances", file_name); comments=true, comment_char='/')

            status.number_items, status.bin_capacity, status.best_solution = data[1:3]
            status.weight = sort(data[4:end], rev=true)
            status.total_accumulated_weight = sum(status.weight)

            # Sort_Descending_Weights
            status.ordered_weight = sortperm(status.weight, rev=true)
            status.L2, status.n_ = LowerBound(status)
            status.permutation = status.ordered_weight[(status.n_+1):end]

            GGA_CGT(status)

            println("Generation:\t", status.generation)
            println("        L2:\t", status.best_solution)
            println("Best sol. :\t", length(status.global_best_solution.bins))
            println("Total time:\t", status.TotalTime)
            println("--------------------------------------------------")


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
        status.global_best_solution = Find_Best_Solution(status)
        
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

        if !is_feasible(individual, status)
            @error "Error found in Generate_Initial_Population"
            exit()
            return true
        end
        
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
    
    update_solution(individual, status)
    
    if !is_feasible(individual, status)
        @error "Error found in FF_n_"
        exit()
    end
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
      
        if stop_criteria_is_met(child2, status)
            return true
        end
        
        push!(children, child2)
        i-=1; j+=2
    end


    # ==========================================================================
    #               Controlled replacement for crossover
    # ==========================================================================
    k = round(Int, p_c/2*P_size - 1)
    status.population[random_individuals[1:k]] = children[1:k]


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
    # aqui para abajo falta (parece que funciona)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # ==========================================================================
    #                   Controlled selection for mutation 
    # ==========================================================================
    sort!(status.population, lt = (a, b) -> a.fitness < b.fitness)

    j = 1
    for i = reverse(1:round(Int,P_size - (p_m*P_size)))
  
        child = Solution(deepcopy(status.population[i].bins), status)

        if !is_feasible(child, status)
            @show i, j
            @error "Error found in Generation 00"
            display(child.bins)
            exit()
        end
        
        if j < P_size*B_size && status.generation - status.population[i].generation < status.life_span
      
            # ==================================================================
            #           Controlled replacement for mutation
            # ==================================================================
            Adaptive_Mutation_RP(child, status.k_cs, true, status)

            if !is_feasible(child, status)
                @show i, j
                @error "Error found in Generation 2"
                exit()
            end

            status.population[j] = child
          
            if stop_criteria_is_met(status.population[j], status)
                return true
            end
          
            j+=1;
        else
            Adaptive_Mutation_RP(child, status.k_ncs, false, status);

            status.population[i] = child

            if !is_feasible(child, status)
                @show i, j
                @error "Error found in Generation 3"
                exit()
            end
            
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

    s =  Solution(bins, status)

    if !is_feasible(s, status)
        @error "Error desde el Gene_Level_Crossover_FFD"
        exit()
    end

    s


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

    if !is_feasible(individual, status)
        @show is_cloned
        @error "Error found in Adaptive_Mutation_RP 1"
        exit()
    end

    pow(a,b) = a^b

    ordered_BinFullness = collect(1:individual.number_of_bins)
   
    if is_cloned
        shuffle!(individual.bins)
    else
        sort!(individual.bins; lt = (a, b) -> a.Bin_Fullness < b.Bin_Fullness)
    end

    # Sort_Ascending_BinFullness

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
 

    RP(individual, findall(free_items), status)

    update_solution(individual, status)
    if !is_feasible(individual, status)
        @error "Error found in Adaptive_Mutation_RP"
        exit()
    end

end



function swap!(bin, weight, id1, id2, free_items, new_free_items, bin_capacity)
    number_free_items = length(free_items)
    p, s = bin.w[id1], bin.w[id2]

    sw = weight[p] + weight[s]

    for i = 1:number_free_items-1
        a = free_items[i]
        if weight[a] >= sw && bin.Bin_Fullness + weight[a]  - sw <= bin_capacity
            bin.Bin_Fullness += weight[a]  - sw
            bin.w[id1] = a

            push!(new_free_items, p, s)
            deleteat!(bin.w, id2)
            deleteat!(free_items, i)
            return true
        end
        

        for j = i+1:number_free_items
            b=free_items[j]
            
            if weight[b] >= sw && bin.Bin_Fullness + weight[b]  - sw <= bin_capacity
                bin.Bin_Fullness += weight[b]  - sw
                bin.w[id1] = b

                push!(new_free_items, p, s)
                deleteat!(bin.w, id2)
                deleteat!(free_items, j)
                return true
            end
            if weight[a] + weight[b] >= sw && bin.Bin_Fullness + weight[b] + weight[a]  - sw <= bin_capacity
                bin.Bin_Fullness += weight[a] + weight[b]   - sw
                
                bin.w[id1] = a
                bin.w[id2] = b

                push!(new_free_items, p, s)
                deleteat!(free_items, (i,j))
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
        for i = 1:length(bin.w)-1
            for j = i+1:length(bin.w)
                update_sol_flag = swap!(bin, weight, i, j, F, new_free_items, status.bin_capacity)
                update_sol_flag && break
            end
            update_sol_flag && break
        end
        
    end

    push!(new_free_items, F...)

    shuffle!(new_free_items)


    for item = new_free_items
        FF(item, individual, status);
    end
    
end

function is_feasible(individual, status)
    used_items = zeros(Int, status.number_items)
    for bin = individual.bins
        used_items[bin.w] .+= 1
    end

    if sum(used_items) != status.number_items || findfirst(a -> a != 1, used_items) != nothing
        @error "Unexpected solution found:"
        @show individual.bins
        @show used_items
        @show status.generation
        return false

    end

    return true
    
end


main()