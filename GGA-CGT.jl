"""
/************************************************************************************************************************
GGA-CGT                                                                                                       *
GROUPING GENETIC ALGORITHM WITH CONTROLLED GENE TRANSMISSION FOR THE BIN PACKING PROBLEM                      *
************************************************************************************************************************/
/************************************************************************************************************************
 Author:  Marcela Quiroz-Castellanos                                                          *
      qc.marcela@gmail.com                                                                            *
         Tecnológico Nacional de México                                                                 *
         Instituto Tecnológico de Ciudad Madero                                                         *
         División de Estudios de Posgrado e Investigación                                               *
         Depto. de Sistemas y Computación                                                               *
************************************************************************************************************************/
/************************************************************************************************************************
 The program excecutes GGA-CGT over a set instances using different configurations                          *
 given by the user. Each configuration represents an independent execution of the GA.                       *
                                                                                *
 Reference:                                                                                                   *
    Quiroz-Castellanos, M., Cruz-Reyes, L., Torres-Jimenez, J.,                                     *
    Gómez, C., Huacuja, H. J. F., & Alvim, A. C. (2015).                                          *
      A grouping genetic algorithm with controlled gene transmission for                                *
      the bin packing problem. Computers & Operations Research, 55, 52-64.                                *
                                                                                                          *
 Input:                                                                                                       *
  File "instances.txt" including the name of the BPP instances to be solved;                            *
  Files including the standard instances to be solve;                                           *
  File "configurations.txt" including the parameter values for each experiment;                           *
                                                                                                              *
 Output:                                                                            *
  A set of files "GGA-CGT_(i).txt" including the experimental results for each                            *
  configuration i in the input, stored in directory: Solutions_GGA-CGT;                               *
   If(save_bestSolution = 1) a set of files HGGA_S_(i)_instance.txt including the                         *
   obtained solution for each instance, for each configuration i, stored in directory: Details_GGA-CGT;           *
************************************************************************************************************************/
"""


# CONSTANTS DEFINING THE SIZE OF THE PROBLEM
const ATTRIBUTES = 5000
const P_size_MAX = 500

# char
#   file[300],
#   nameC[300];

# int
#   is_optimal_solution,
#    save_bestSolution,
#   generation,
#    repeated_fitness,
#    max_gen,
#    life_span,
#    P_size,
#    seed;

# long int
#   i,
#    j,
#    conf,
#   number_items,
#    bin_capacity,
#   best_solution,
#    n_,
#    L2,
#    bin_i,
#    higher_weight,
#    lighter_weight,
#   ordered_weight[ATTRIBUTES],
#    weight[ATTRIBUTES],
#    permutation[ATTRIBUTES],
#    items_auxiliary[ATTRIBUTES],
#    ordered_population[P_size_MAX],
#    best_individuals[P_size_MAX],
#    random_individuals[P_size_MAX];

# float
#   p_m,
#   p_c,
#   k_ncs,
#   k_cs,
#   B_size,
#   TotalTime;

# long double
#   total_accumulated_weight,
#    weight1[ATTRIBUTES],
#   _p_;

# clock_t
#   start,
#   end;

# FILE  *output,
#     *input_Configurations,
#       *input_Instances;

struct SOLUTION
  linked_list L;
  double Bin_Fullness;
end

# SOLUTION  global_best_solution[ATTRIBUTES],
#       population[P_size_MAX][ATTRIBUTES],
#          children[P_size_MAX][ATTRIBUTES];

# Initial seeds for the random number generation
# int   seed_emptybin,
#     seed_permutation;

# GA COMPONENTS
# long int Generate_Initial_Population();
# long int Generation();
# void  Gene_Level_Crossover_FFD(long int, long int, long int);
# void  Adaptive_Mutation_RP(long int, float, int);
# void  FF_n_(int); # First Fit with ñ pre-allocated-items (FF-ñ)
# void  RP(long int, long int &, long int[], long int); # Rearrangement by Pairs

# # BPP Procedures
# void  FF(long int, SOLUTION[], long int&, long int, int);
# void  LowerBound();

# # Auxiliary functions
# void  Find_Best_Solution();
# void  Sort_Ascending_IndividualsFitness();
# void  Sort_Descending_Weights(long int[], long int);
# void  Sort_Ascending_BinFullness(long int[], long int);
# void  Sort_Descending_BinFullness(long int [], long int);
# void  Sort_Random(long int[], long int, int);
# void  Copy_Solution(SOLUTION[], SOLUTION[], int);
# void  Clean_population();
# long int  Used_Items(long int, long int, long int[]);
# void  Adjust_Solution(long int);
# long int  LoadData();
# void  WriteOutput();
# void  sendtofile(SOLUTION[]);

# # Pseudo-random number generator functions
# int get_rand_ij(int *, int, int);
# int   get_rand(int *, int);
# float randp(int *);
# int trand();

import CSV

function main ()

    # char  aux[300], nombreC[300], mystring[300];
    try
        mkdir("Solutions_GGA-CGT");
        mkdir("Details_GGA-CGT");
    end

    # READING EACH CONFIGURATION IN FILE "configurations.txt", CONTAINING THE PARAMETER VALUES FOR EACH EXPERIMENT
    input_Configurations = nothing
    try
        input_Configurations = CSV.File("configurations.txt")
    catch
        print("\n INVALID FILE");
        return false
    end
     
    for nconf = length(input_Configurations)
        conf = input_Configurations.conf[1]
        P_size = input_Configurations.P_size[1]
        max_gen = input_Configurations.max_gen[1]
        p_m = input_Configurations.p_m[1]
        p_c = input_Configurations.p_c[1]
        k_ncs = input_Configurations.k_ncs[1]
        k_cs = input_Configurations.k_cs[1]
        B_size = input_Configurations.B_size[1]
        life_span = input_Configurations.life_span[1]
        seed = input_Configurations.seed[1]
        save_bestSolution = input_Configurations.save_bestSolution[1]
        
        nameC = string("Solutions_GGA-CGT/GGA-CGT_(", conf, ").txt");
        open(nameC,"w+") do output
            write(output , "CONF\t|P|\tmax_gen\tn_m\tn_c\tk1(non-cloned_solutions)\tk2(cloned_solutions)\t|B|\tlife_span\tseed");
            write(output, "\n$conf\t$P_size\t$max_gen\t$p_m\tp_c\t$k_ncs\t$k_cs\t$B_size\t$life_span\t$seed");
            write(output, "\nInstancias \t L2 \t Bins \t Gen \t Time");
        end
        
        # input_Instances = nothing
        # # READING FILE "instances.txt" CONTAINING THE NAME OF BPP INSTANCES TO PROCESS
        # try
        #     input_Instances = open("instances.txt","rt")
        # catch
        #     print("\n INVALID FILE");
        #     # getch();
        #     return 1;
        # end

        for file in readlines("instances.txt")
            LoadData();

            for i = 1:number_items
                ordered_weight[i] = i;
            end

            Sort_Descending_Weights(ordered_weight, number_items);
            LowerBound();
            seed_permutation = seed;
            seed_emptybin = seed;
            
            for i = 1:P_size
                ordered_population[i] = i;
                random_individuals[i] = i;
                best_individuals[i] = i;
            end

            Clean_population();
            is_optimal_solution = 0;
            generation = 0;
            j = n_
            i = 0

            for j = 1:number_items 
                permutation[i] = ordered_weight[j];
                i+=1
            end
             
            repeated_fitness = 0;
            # procedure GGA-CGT
            start = clock();
            if !Generate_Initial_Population()  # Generate_Initial_Population() returns 1 if an optimal solution was found
                for generation = 1:max_gen)
                    if Generation() # Generation() returns 1 if an optimal solution was found
                        break;
                    end
                    Find_Best_Solution();
                end
                if(!is_optimal_solution) # is_optimal_solution is 1 if an optimal solution was printed before
                    end_time = clock();
                    TotalTime = (end_time - start) / (CLOCKS_PER_SEC * 1.0);
                    Find_Best_Solution();
                    WriteOutput();
                end
            end

        end
        fclose(input_Instances);
    end
    fclose(input_Configurations);
    print("\n\n\tEnd of process");
    getchar();

  return 0;
end



"""
/************************************************************************************************************************
 To generate an initial population P of individuals with FF-ñ packing heuristic.                            *
  population[i][number_items].Bin_Fullness: Saves the fitness of the solution i                                         *
  population[i][number_items + 1].Bin_Fullness:                                                                         *
  Saves the total number of bins in the solution i                                                                      *
  population[i][number_items + 2].Bin_Fullness:                                                                         *
  Saves the generation in which the solution i was generated                                                          *
  population[i][number_items + 3].Bin_Fullness:                                                                         *
  Saves the number of bins in the solution i that are fully at 100%                                                     *
  population[i][number_items + 4].Bin_Fullness:                                                   *
  Saves the fullness of the bin with the highest avaliable capacity in the solution i                                 *
 Output:                                                                                                        *
  (1) when it finds a solution for which the size matches the L2 lower bound                              *
   (0) otherwise                                                                      *
************************************************************************************************************************/
"""
function Generate_Initial_Population()
    for i = 1:P_size
        FF_n_(i);
        population[i][number_items + 2].Bin_Fullness = generation;
        population[i][number_items].Bin_Fullness /= population[i][number_items + 1].Bin_Fullness;
        
        if(population[i][number_items + 1].Bin_Fullness == L2)
            end_time = clock();
            Copy_Solution(global_best_solution, population[i], 0);
            global_best_solution[number_items].Bin_Fullness = population[i][number_items].Bin_Fullness;;
            global_best_solution[number_items + 2].Bin_Fullness = generation;
            global_best_solution[number_items + 1].Bin_Fullness = population[i][number_items + 1].Bin_Fullness;
            global_best_solution[number_items + 3].Bin_Fullness = population[i][number_items + 3].Bin_Fullness;
            TotalTime = (end_time - start) / (CLOCKS_PER_SEC * 1.0);
            WriteOutput();
            is_optimal_solution = 1;
            return 1;
        end
    
    end
    
    return 0;
end



"""
/************************************************************************************************************************
 To apply the reproduction technique: Controlled selection and Controlled replacement.                          *
 Output:                                                                                                        *
  (1) when it finds a solution for which the size matches the L2 lower bound                              *
   (2) if more than 0.1*P_size individuals (solutions) have duplicated-fitness                            *
   (0) otherwise                                                                      *
************************************************************************************************************************/
"""
function Generation()

   # long int
   #    f1,
   #    f2,
   #       h,
    #   k;

   """
   /*****************************************************************************************************-
   ********************************-Controlled selection for crossover************************************
   ****************************************************************************************************-*/
   """
    Sort_Ascending_IndividualsFitness();
    # if(generation > 1 && repeated_fitness > 0.1*P_size)
    #   return(2);
    Sort_Random(random_individuals,0,(int)(P_size-(int)(P_size*B_size)));
    Sort_Random(best_individuals,(1-p_c)*P_size,P_size);
    k = 0;
    h = P_size - 1;
    i = P_size - 1; j = 0
   
    while (i > P_size - (p_c/2*P_size))
        f1 = best_individuals[h-=1];
        f2 = random_individuals[k+=1];
      
        if f2 == f1
            f1 = best_individuals[h-=1];
        end
      
        Gene_Level_Crossover_FFD(ordered_population[f1], ordered_population[f2], j);
        children[j][number_items + 2].Bin_Fullness = generation + 1;
        children[j][number_items].Bin_Fullness /= children[j][number_items+1].Bin_Fullness;
      
        if(children[j][number_items + 1].Bin_Fullness == L2)  
            end_time = clock();
            Copy_Solution(global_best_solution, children[j], 0);
            global_best_solution[number_items].Bin_Fullness = children[j][number_items].Bin_Fullness;;
            global_best_solution[number_items + 2].Bin_Fullness = generation+1;
            global_best_solution[number_items + 1].Bin_Fullness = children[j][number_items + 1].Bin_Fullness;
            global_best_solution[number_items + 3].Bin_Fullness = children[j][number_items + 3].Bin_Fullness;
            TotalTime = (end_time - start) / (CLOCKS_PER_SEC * 1.0);
            WriteOutput();
            is_optimal_solution = 1;
            return 1;
        end
      
        Gene_Level_Crossover_FFD(ordered_population[f2], ordered_population[f1], j+1);
        children[j+1][number_items + 2].Bin_Fullness = generation + 1;
        children[j+1][number_items].Bin_Fullness /= children[j+1][number_items+1].Bin_Fullness;
      
        if(children[j+1][number_items + 1].Bin_Fullness == L2)   
            end_time = clock();
            Copy_Solution(global_best_solution, children[j+1], 0);
            global_best_solution[number_items].Bin_Fullness = children[j+1][number_items].Bin_Fullness;;
            global_best_solution[number_items + 2].Bin_Fullness = generation+1;
            global_best_solution[number_items + 1].Bin_Fullness = children[j+1][number_items + 1].Bin_Fullness;
            global_best_solution[number_items + 3].Bin_Fullness = children[j+1][number_items + 3].Bin_Fullness;
            TotalTime = (end_time - start) / (CLOCKS_PER_SEC * 1.0);
            WriteOutput();
            is_optimal_solution = 1;
            return 1;
        end
        
        i-=1; j+=2
    end

    """
    /*****************************************************************************************************-
    ********************************-Controlled replacement for crossover**********************************
    ****************************************************************************************************-*/
    """
    k = 0;
    for j = 1:p_c/2*P_size - 1
        Copy_Solution(population[ordered_population[random_individuals[k]]], children[j], 1);
        k+=1
    end

    k = 0;

   for(i = P_size - 1; i > P_size - (p_c/2*P_size); i-=1, j+=1)
        while population[ordered_population[k]][number_items + 2].Bin_Fullness == generation + 1
            k+=1;
        end
        Copy_Solution(population[ordered_population[k+=1]], children[j], 1);
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
  
        if(j < P_size*B_size && generation+ 1 - population[ordered_population[i]][number_items + 2].Bin_Fullness < life_span)
      
            """
            /*****************************************************************************************************-
            **********************************Controlled replacement for mutation**********************************
            ****************************************************************************************************-*/
            """
            Copy_Solution(population[ordered_population[j]], population[ordered_population[i]], 0);
            Adaptive_Mutation_RP(ordered_population[j], k_cs, 1);
            population[ordered_population[j]][number_items + 2].Bin_Fullness = generation + 1;
            population[ordered_population[j]][number_items].Bin_Fullness /= population[ordered_population[j]][number_items+1].Bin_Fullness;
          
            if(population[ordered_population[j]][number_items + 1].Bin_Fullness == L2)
                end_time = clock();
                Copy_Solution(global_best_solution, population[ordered_population[j]], 0);
                global_best_solution[number_items].Bin_Fullness = population[ordered_population[j]][number_items].Bin_Fullness;;
                global_best_solution[number_items + 2].Bin_Fullness = generation+1;
                global_best_solution[number_items + 1].Bin_Fullness = population[ordered_population[j]][number_items + 1].Bin_Fullness;
                global_best_solution[number_items + 3].Bin_Fullness = population[ordered_population[j]][number_items + 3].Bin_Fullness;
                TotalTime = (end_time - start) / (CLOCKS_PER_SEC * 1.0);
                WriteOutput();
                is_optimal_solution = 1;
                return 1;
            end
          
            j+=1;
        else
            Adaptive_Mutation_RP(ordered_population[i], k_ncs, 0);
            population[ordered_population[i]][number_items + 2].Bin_Fullness = generation + 1;
            population[ordered_population[i]][number_items].Bin_Fullness /= population[ordered_population[i]][number_items+1].Bin_Fullness;
            if population[ordered_population[i]][number_items+1].Bin_Fullness == L2
                end_time = clock();
                Copy_Solution(global_best_solution, population[ordered_population[i]], 0);
                global_best_solution[number_items].Bin_Fullness = population[ordered_population[i]][number_items].Bin_Fullness;;
                global_best_solution[number_items + 2].Bin_Fullness = generation+1;
                global_best_solution[number_items + 1].Bin_Fullness = population[ordered_population[i]][number_items + 1].Bin_Fullness;
                global_best_solution[number_items + 3].Bin_Fullness = population[ordered_population[i]][number_items + 3].Bin_Fullness;
                TotalTime = (end_time - start) / (CLOCKS_PER_SEC * 1.0);
                WriteOutput();
                is_optimal_solution = 1;
                return 1;
            end
        end
    end
    
    return 0;
end



"""
/************************************************************************************************************************
 To recombine two parent solutions producing a child solution.                                          *
 Input:                                                                                                         *
  The positions in the population of the two parent solutions: father_1 and father_2                        *
  The position in the set of children of the child solution: child                                    *
************************************************************************************************************************/
"""
function Gene_Level_Crossover_FFD(long int father_1, long int father_2, long int child)

  # long int  k,
 #      counter,
 #      k2 = 0,
 #      ban = 1,
 #         items[ATTRIBUTES] = {0},
 #         free_items[ATTRIBUTES] = {0};
 #   children[child][number_items + 4].Bin_Fullness = bin_capacity;

    if(population[father_1][number_items + 1].Bin_Fullness > population[father_2][number_items + 1].Bin_Fullness)
        counter = population[father_1][number_items + 1].Bin_Fullness;
    else
        counter = population[father_2][number_items + 1].Bin_Fullness;
    end

    long int *random_order1 = new long int [counter];
    long int *random_order2 = new long int [counter];


    for(k = 1:counter)
        random_order1[k] = k;
        random_order2[k] = k;
    end

    Sort_Random(random_order1,0, population[father_1][number_items + 1].Bin_Fullness);
    Sort_Random(random_order2,0, population[father_2][number_items + 1].Bin_Fullness);
    Sort_Descending_BinFullness(random_order1, father_1);
    Sort_Descending_BinFullness(random_order2, father_2);

    for(k = 1:population[father_1][number_items + 1].Bin_Fullness)
   
        if(population[father_1][random_order1[k]].Bin_Fullness >= population[father_2][random_order2[k]].Bin_Fullness)
            
            ban = Used_Items(father_1, random_order1[k], items);
            if (ban == 1)
                children[child][k2].L.clone_linked_list(population[father_1][random_order1[k]].L);
                children[child][k2+=1].Bin_Fullness = population[father_1][random_order1[k]].Bin_Fullness;
            
                if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
                    children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
                end
            end
            
            if(population[father_2][random_order2[k]].Bin_Fullness > 0)
                ban = Used_Items(father_2, random_order2[k], items);
                if (ban == 1)
                    children[child][k2].L.clone_linked_list(population[father_2][random_order2[k]].L);
                    children[child][k2+=1].Bin_Fullness = population[father_2][random_order2[k]].Bin_Fullness;

                    if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
                        children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
                    end
                end
            end
        else
            if(population[father_2][random_order2[k]].Bin_Fullness > 0)
                ban = Used_Items(father_2, random_order2[k], items);
                if (ban == 1)
                    children[child][k2].L.clone_linked_list(population[father_2][random_order2[k]].L);
                    children[child][k2+=1].Bin_Fullness = population[father_2][random_order2[k]].Bin_Fullness;
                    
                    if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
                        children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
                    end
                end
            end

            ban = Used_Items(father_1, random_order1[k], items);
            if (ban == 1)
                children[child][k2].L.clone_linked_list(population[father_1][random_order1[k]].L);
                children[child][k2+=1].Bin_Fullness = population[father_1][random_order1[k]].Bin_Fullness;
                if(children[child][k2-1].Bin_Fullness < children[child][number_items + 4].Bin_Fullness)
                    children[child][number_items + 4].Bin_Fullness = children[child][k2-1].Bin_Fullness;
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
            children[child][number_items].Bin_Fullness += pow((children[child][k].Bin_Fullness / bin_capacity), 2);
        end
    end
    children[child][number_items+1].Bin_Fullness = k2;
    free(random_order1);
    free(random_order2);
end



"""
/************************************************************************************************************************
 To produce a small modification in a solution.                                                 *
 Input:                                                                                                         *
  The position in the population of the solution to mutate: individual                                *
  The rate of change to calculate the number of bins to eliminate: k                                  *
  A value that indicates if the solution was cloned: is_cloned                                      *
************************************************************************************************************************/
"""
function Adaptive_Mutation_RP(long int individual, float k, int is_cloned)

  # long int
  #        number_bins,
  #        i,
  #        lightest_bin = 0,
  #        number_free_items = 0,
  #        free_items[ATTRIBUTES] = {0},
  #        ordered_BinFullness[ATTRIBUTES] = {0};
  #  node *p;

    for(i = 1:population[individual][number_items + 1].Bin_Fullness)
        ordered_BinFullness[i] = i;
    end
   
    if(is_cloned)
        Sort_Random(ordered_BinFullness,0, population[individual][number_items + 1].Bin_Fullness);
    end

    Sort_Ascending_BinFullness(ordered_BinFullness, individual);
    i = 1;
    
    while(population[individual][ordered_BinFullness[i]].Bin_Fullness < bin_capacity && i < population[individual][number_items + 1].Bin_Fullness)
        i+=1;
    end    
    
    _p_ = 1 / (float)(k);
    number_bins = (long int)ceil(i*((2 - i/population[individual][number_items + 1].Bin_Fullness) / pow(i,_p_))*(1 - ((double)get_rand(&seed_emptybin,(long int)ceil((1/pow(i,_p_))*100))/100)));
    
    for i = 1:number_bins
    
        p = population[individual][ordered_BinFullness[lightest_bin]].L.first;

        while(p != NULL)
            free_items[number_free_items+=1] = p->data;
            p = p->next;
        end

        population[individual][ordered_BinFullness[lightest_bin]].L.free_linked_list();
        population[individual][ordered_BinFullness[lightest_bin]].Bin_Fullness = 0;
        lightest_bin+=1;
    end
 
    population[individual][number_items + 1].Bin_Fullness -= number_bins;
    number_bins = population[individual][number_items + 1].Bin_Fullness;
 
    Adjust_Solution(individual);
    RP(individual, number_bins, free_items, number_free_items);
 
    population[individual][number_items+1].Bin_Fullness = number_bins;
end



"""
/************************************************************************************************************************
 To generate a random BPP solution with the ñ large items packed in separate bins.                          *
 Input:                                                                                                         *
  The position in the population of the new solution: individual                                    *
************************************************************************************************************************/
"""
function FF_n_(int individual)
  # long int
  #     i,
  #     j = 0,
  #     total_bins = 0;

    bin_i = 0;
    population[individual][number_items + 3].Bin_Fullness = 0;
    population[individual][number_items + 4].Bin_Fullness = bin_capacity;
    if (n_ > 0)
        for i = 1:n_
            population[individual][i].Bin_Fullness = weight[ordered_weight[i]];
            population[individual][i].L.insert(ordered_weight[i]);
            total_bins+=1;

            if(population[individual][i].Bin_Fullness < population[individual][number_items + 4].Bin_Fullness)
                population[individual][number_items + 4].Bin_Fullness = population[individual][i].Bin_Fullness;
            end
        end
    
        i = number_items - i;
        Sort_Random(permutation,0, i);
        for(j = 0; j < i-1; j+=1)
            FF(permutation[j], population[individual], total_bins, bin_i,0);
        end
        FF(permutation[j], population[individual], total_bins, bin_i,1);
   
    else
        Sort_Random(permutation,0, number_items);
        for(j = 0; j < number_items-1; j+=1)
            FF(permutation[j], population[individual], total_bins, bin_i,0);
            FF(permutation[j], population[individual], total_bins, bin_i,1);
        end
   end

  population[individual][number_items + 1].Bin_Fullness = total_bins;
end



"""
/************************************************************************************************************************
 To reinsert free items into an incomplete BPP solution.                                            *
 Input:                                                                                                                 *
  The position in the population of the incomplete solution where the free items must be reinserted: individual     *
   The number of bins of the partial_solution: b                                                *
   A set of free items to be reinserted into the partial_solution: F                                    *
   The number of free items of F: number_free_items                                             *
************************************************************************************************************************/
"""
function RP(long int individual, long int &b, long int F[], long int number_free_items)

  # long int
  #     i,
  #        k,
  #        k2,
  #        ban,
  #        sum = 0,
  #        total_free = 0,
  #        ordered_BinFullness[ATTRIBUTES] = {0},
  #     *new_free_items = new long int [2];

   # node   *p,
   #    *s,
   #       *aux;

   higher_weight = weight[F[0]];
   lighter_weight = weight[F[0]];
   bin_i = b;
   population[individual][number_items].Bin_Fullness = 0;
   population[individual][number_items + 3].Bin_Fullness = 0;
   population[individual][number_items + 4].Bin_Fullness = bin_capacity;

    for (i = 0; i < b; i+=1 )
        ordered_BinFullness[i] = i;
    end
    
    Sort_Random(ordered_BinFullness,0, b);
    Sort_Random(F,0, number_free_items);

    for (i = 0; i < b; i+=1)
        sum = (long int)population[individual][ordered_BinFullness[i]].Bin_Fullness;

        p = population[individual][ordered_BinFullness[i]].L.first;
        while(p->next != NULL)
            ban = 0;
            aux = p;
            s = p->next;
            while(s != NULL)       
                for (k = 0; k < number_free_items - 1; k+=1)
                    if(i == b-1)
                        if(weight[F[k]] > higher_weight)
                            higher_weight = weight[F[k]];
                        end
                    end
                    for (k2 = k + 1; k2 < number_free_items; k2+=1)
                        if(weight[F[k]] >= weight[p->data] + weight[s->data] && ((sum - (weight[p->data] + weight[s->data]) + (weight[F[k]]) <= bin_capacity)))
                            sum = sum - (weight[p->data] + weight[s->data]) + (weight[F[k]]);
                            new_free_items[0] = p->data;
                            new_free_items[1] = s->data;
                            p->data = F[k];
                            aux->next = s->next;
                            free(s);
                            
                            if(population[individual][ordered_BinFullness[i]].L.last == s)
                                population[individual][ordered_BinFullness[i]].L.last = aux;
                            end
                            
                            population[individual][ordered_BinFullness[i]].L.num-=1;
                            F[k] = new_free_items[0];
                            F[number_free_items + total_free] = new_free_items[1];
                            total_free+=1;
                            ban = 1;
                            break;
                        end

                        if(weight[F[k2]] >= weight[p->data] + weight[s->data] && ((sum - (weight[p->data] + weight[s->data]) + (weight[F[k2]]) <= bin_capacity)))
                            sum = sum - (weight[p->data] + weight[s->data]) + (weight[F[k2]]);
                            new_free_items[0] = p->data;
                            new_free_items[1] = s->data;
                            p->data = F[k2];
                            aux->next = s->next;
                            free(s);
                            
                            if(population[individual][ordered_BinFullness[i]].L.last == s)
                                population[individual][ordered_BinFullness[i]].L.last = aux;
                            end

                            population[individual][ordered_BinFullness[i]].L.num-=1;
                            F[k2] = new_free_items[0];
                            F[number_free_items + total_free] = new_free_items[1];
                            total_free+=1;
                            ban = 1;
                            break;

                        end

                        if((weight[F[k]] + weight[F[k2]] > weight[p->data] + weight[s->data]) || ((weight[F[k]] + weight[F[k2]] == weight[p->data] + weight[s->data]) && !(weight[F[k]] == weight[p->data] || weight[F[k]] == weight[s->data])))
                          
                            if(sum - (weight[p->data] + weight[s->data]) + (weight[F[k]] + weight[F[k2]]) > bin_capacity)
                                break;
                            end
                        
                            sum = sum - (weight[p->data] + weight[s->data]) + (weight[F[k]] + weight[F[k2]]);
                            new_free_items[0] = p->data;
                            new_free_items[1] = s->data;
                            p->data = F[k];
                            s->data = F[k2];
                            F[k] = new_free_items[0];
                            F[k2] = new_free_items[1];
                            if(sum == bin_capacity)
                                ban = 1;
                                break;
                            end
                        end
                    end

                    if(ban)
                        break;
                    end

                end

                if(ban)
                    break;
                end

                aux = s;
                s = s->next;

            end
            
            if(ban)
                break;
            end
                p = p->next;
        end
          
        population[individual][ordered_BinFullness[i]].Bin_Fullness = sum;
      
        if(population[individual][ordered_BinFullness[i]].Bin_Fullness < population[individual][number_items + 4].Bin_Fullness)
            population[individual][number_items + 4].Bin_Fullness = population[individual][ordered_BinFullness[i]].Bin_Fullness;
        end
        
        if(population[individual][ordered_BinFullness[i]].Bin_Fullness == bin_capacity)
            population[individual][number_items + 3].Bin_Fullness+=1;
        elseif(population[individual][ordered_BinFullness[i]].Bin_Fullness + weight[ordered_weight[number_items-1]] <= bin_capacity)
            if(ordered_BinFullness[i] < bin_i)
                bin_i = ordered_BinFullness[i];
            end
        end
    end

    for(i = 0; i < bin_i; i+=1)
        population[individual][number_items].Bin_Fullness += pow((population[individual][i].Bin_Fullness / bin_capacity), 2);
    end
   
    free(new_free_items);
    number_free_items += total_free;

    if(higher_weight < .5*bin_capacity)
        Sort_Random(F,0, number_free_items);
    else
        Sort_Descending_Weights(F, number_free_items);
        lighter_weight = weight[F[number_free_items-1]];
    end

    if(lighter_weight > bin_capacity - population[individual][number_items + 4].Bin_Fullness)
        for(i = bin_i; i < b; i+=1)
            population[individual][number_items].Bin_Fullness += pow((population[individual][i].Bin_Fullness / bin_capacity), 2);
        end
        bin_i = b;
    end

    for(i = 0; i < number_free_items-1; i+=1)
        FF(F[i], population[individual], b, bin_i, 0);
    end
    
    FF(F[i], population[individual], b, bin_i, 1);
end



"""
/************************************************************************************************************************
 To insert an item into an incomplete BPP solution.                                             *
 Input:                                                                                                         *
   An item to be inserted into the individual: item                                             *
  An incomplete chromosome where the item must be inserted: individual                                      *
   The number of bins of the individual: total_bins                                             *
   The first bin that could have sufficient available capacity to store the item: beginning                   *
   A value that indicates if it is the last item to be stored into the individual: is_last                    *
************************************************************************************************************************/
"""
function FF(long int item, SOLUTION individual[], long int &total_bins, long int beginning, int is_last)

  # long int  i;

    if(!is_last && weight[item] > (bin_capacity - (long int)individual[number_items + 4].Bin_Fullness))
        i = total_bins;
    else
        for(i = beginning; i < total_bins; i+=1)
            if((long int)individual[i].Bin_Fullness + weight[item] <= bin_capacity)
          
                individual[i].Bin_Fullness += weight[item];
                individual[i].L.insert(item);
                if((long int)individual[i].Bin_Fullness == bin_capacity)
                    individual[number_items + 3].Bin_Fullness+=1;
                end

                if(is_last)
                    for(i=i; i < total_bins; i+=1)
                        individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
                    end
                    return;
                end

                if((long int)individual[i].Bin_Fullness + weight[ordered_weight[number_items-1]] > bin_capacity && i == bin_i)
                    bin_i+=1;
                    individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
                end
                return;
            end
        end
        if(is_last)
            individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
        end
    end
    individual[i].Bin_Fullness += weight[item];
    individual[i].L.insert(item);
    
    if(individual[i].Bin_Fullness < individual[number_items + 4].Bin_Fullness)
        individual[number_items + 4].Bin_Fullness = individual[i].Bin_Fullness;
    end
    
    if(is_last)
      individual[number_items].Bin_Fullness += pow((individual[i].Bin_Fullness / bin_capacity), 2);
    end
    
    total_bins+=1;
end



"""
/************************************************************************************************************************
 To calculate the lower bound L2 of Martello and Toth and the ñ large items n_                            *
************************************************************************************************************************/
"""
function LowerBound()

    long int k, m, i, j, aux1, aux2;
    long double sjx=0, sj2=0, sj3=0;
    long int jx=0, cj12, jp=0, jpp=0, cj2;

    while(weight[ordered_weight[jx]] > bin_capacity/2 && jx < number_items)
        jx+=1;
    end
   
    n_ = jx;
    if(jx == number_items)
        L2 = jx;
        return;
    end
   
    if(jx == 0)
        if(fmod(total_accumulated_weight,bin_capacity) >= 1)
            L2 = (long int)ceil(total_accumulated_weight / bin_capacity);
        else
            L2 = (long int)(total_accumulated_weight / bin_capacity);
            return;
        end
    else
        cj12 = jx;
        for(i=jx; i < number_items; i+=1)
            sjx += weight[ordered_weight[i]];
        end

        jp = jx;
        
        for(i = 0; i < jx; i+=1)
            if(weight[ordered_weight[i]] <= bin_capacity - weight[ordered_weight[jx]])
                jp = i;
                break;
            end
        end

        cj2 = jx - jp;
        for(i=jp; i <= jx-1; i+=1)
            sj2 += weight[ordered_weight[i]];
        end

        jpp = jx;
        sj3 = weight[ordered_weight[jpp]];
        ordered_weight[number_items] = number_items;
        weight[number_items]=0;
        
        while(weight[ordered_weight[jpp+1]]==weight[ordered_weight[jpp]])
            jpp+=1;
            sj3 += weight[ordered_weight[jpp]];
        end
 
        L2 = cj12;

        while true
          
            if(fmod((sj3 + sj2),bin_capacity) >= 1)
                aux1 = (long int)ceil((sj3 + sj2)/bin_capacity - cj2);
            else
                aux1 = (long int)((sj3 + sj2)/bin_capacity - cj2);
            end

            if(L2 < (cj12 + aux1))
                L2 = cj12 + aux1;
            end

            jpp+=1;
            if(jpp < number_items)
                sj3 += weight[ordered_weight[jpp]];
                
                while(weight[ordered_weight[jpp+1]] == weight[ordered_weight[jpp]])
                  jpp+=1;
                  sj3 += weight[ordered_weight[jpp]];
                end

                while(jp > 0 && weight[ordered_weight[jp-1]] <= bin_capacity - weight[ordered_weight[jpp]])
                    jp-=1;
                    cj2+=1;
                    sj2 += weight[ordered_weight[jp]];
                end
            
            end
        
            if(fmod((sjx + sj2),bin_capacity) >= 1)
                aux2 = (long int)ceil((sjx + sj2) / bin_capacity - cj2 );
            else
                aux2 = (long int)((sjx + sj2) / bin_capacity - cj2 );
            end

            if(jpp <= number_items || (cj12 + aux2) > L2)
                break
            end
        end  
    end
end



"""
/************************************************************************************************************************
 To find the solution with the highest fitness of the population and update the global_best_solution              *
************************************************************************************************************************/
"""
function Find_Best_Solution()

  # long int  i,
  #     best_individual = 0;
    for(i = 0; i < P_size; i+=1)
        if(population[i][number_items].Bin_Fullness > population[best_individual][number_items].Bin_Fullness)
            best_individual = i;
        end
    end

    if(generation + 1 > 1)
        if(population[best_individual][number_items].Bin_Fullness > global_best_solution[number_items].Bin_Fullness)
            Copy_Solution(global_best_solution, population[best_individual], 0);
        end
    else
        Copy_Solution(global_best_solution, population[best_individual], 0);
    end
end



"""
/************************************************************************************************************************
 To sort the individuals of the population in ascending order of their fitness                              *
************************************************************************************************************************/
"""
function Sort_Ascending_IndividualsFitness()

  # long int i,
  #      k = P_size - 1,
  #      i2 = 0,
  #    aux,
  #      ban = 1;

    while(ban)
        ban = 0;
        for(i = i2; i < k; i+=1)
            if(population[ordered_population[i]][number_items].Bin_Fullness > population[ordered_population[i+1]][number_items].Bin_Fullness)
                aux = ordered_population[i];
                ordered_population[i] = ordered_population[i+1];
                ordered_population[i+1] = aux;
                ban = 1;
            else if(population[ordered_population[i]][number_items].Bin_Fullness == population[ordered_population[i+1]][number_items].Bin_Fullness)
                aux = ordered_population[i+1];
                ordered_population[i+1] = ordered_population[i2];
                ordered_population[i2] = aux;
                i2+=1;
            end
        end
        k-=1;
    end
    repeated_fitness = i2;
end



"""
/************************************************************************************************************************
 To sort the bins of a solution in ascending order of their filling                                   *
 Input:                                                                                                                 *
  An array to save the order of the bins: ordered_BinFullness                                       *                                                                                                                 *
  The position in the population of the solution: individual                                        *
************************************************************************************************************************/
"""
function Sort_Ascending_BinFullness(long int ordered_BinFullness[], long int individual)
  # long int  m,
  #     k,
  #        temporary_variable,
  #        ban = 1;

    k = population[individual][number_items + 1].Bin_Fullness - 1;
    while(ban)
        ban = 0;
        for(m = 0; m < k; m+=1)
            if(population[individual][ordered_BinFullness[m]].Bin_Fullness > population[individual][ordered_BinFullness[m+1]].Bin_Fullness)
                temporary_variable = ordered_BinFullness[m];
                ordered_BinFullness[m] = ordered_BinFullness[m + 1];
                ordered_BinFullness[m + 1] = temporary_variable;
                ban = 1;
            end
        end
        k-=1;
    end
end



"""
/************************************************************************************************************************
 To sort the bins of a solution in descending order of their filling                                    *
 Input:                                                                                                                 *
  An array to save the order of the bins: ordered_BinFullness                                       *                                                                                                                 *
  The position in the population of the solution: individual                                        *
************************************************************************************************************************/
"""
function Sort_Descending_BinFullness(long int ordered_BinFullness[], long int individual)

  # long int  m,
  #     k,
  #        temporary_variable,
  #        ban = 1;

    k = population[individual][number_items + 1].Bin_Fullness - 1;
    while(ban)      
        ban = 0;
        for(m = 0; m < k; m+=1)
            if(population[individual][ordered_BinFullness[m]].Bin_Fullness < population[individual][ordered_BinFullness[m+1]].Bin_Fullness)          
                temporary_variable = ordered_BinFullness[m];
                ordered_BinFullness[m] = ordered_BinFullness[m + 1];
                ordered_BinFullness[m + 1] = temporary_variable;
                ban = 1;
            end
        end
        k-=1;
    end
end



"""
/************************************************************************************************************************
 To sort the elements between the positions [k] and [n] of an array in random order                         *
 Input:                                                                                                                 *
  The array to be randomized: random_array                                                    *                                                                                                                 *
  The initial random position: k                                                          *
    The final random position: n                                                            *
************************************************************************************************************************/
"""
function Sort_Random(long int random_array[], long int k, int n)

   # long int i,
   #      aux,
   #      random_number;

    for(i = n - 1; i >= k; i-=1)
        random_number = k + get_rand(&seed_permutation,n-k) - 1;
        aux = random_array[random_number];
        random_array[random_number] = random_array[i];
        random_array[i] = aux;
        if(weight[random_array[i]] < lighter_weight)
            lighter_weight = weight[random_array[i]];
        end
        
        if(weight[random_array[random_number]] < lighter_weight)
            lighter_weight = weight[random_array[random_number]];
        end
    end
end



"""
/************************************************************************************************************************
 To sort a set of items in descending order of their weights                                        *
 Input:                                                                                                                 *
  An array to save the order of the items: ordered_weight                                         *                                                                                                                 *
  The number of items in the set: n                                                       *
************************************************************************************************************************/
"""
function Sort_Descending_Weights(long int ordered_weight[], long int n)
  # long int  m,
  #     k,
  #        temporary_variable,
  #        ban = 1;

    k = n - 1;
    while(ban)
        ban = 0;
        for(m = 0; m < k; m+=1)
            if(weight[ordered_weight[m]] < weight[ordered_weight[m+1]])
                temporary_variable = ordered_weight[m];
                ordered_weight[m] = ordered_weight[m + 1];
                ordered_weight[m + 1] = temporary_variable;
                ban = 1;
            end
        end
        k-=1;
    end
end



"""
/************************************************************************************************************************
 To copy solution2 in solution                                                            *
 Input:                                                                                                                 *
  A solution to save the copied solution: solution                                              *                                                                                                                 *
  The solution to be copied: solution2                                                        *
  A value that indicates if the copied solution must be deleted: delete_solution2                         *
************************************************************************************************************************/
"""
function Copy_Solution(SOLUTION solution[], SOLUTION solution2[], int delete_solution2)
 
  # long int  j;

   for(j = 0; j < solution2[number_items + 1].Bin_Fullness; j+=1)
   
        solution[j].Bin_Fullness = solution2[j].Bin_Fullness;
        solution[j].L.clone_linked_list(solution2[j].L);
        if(delete_solution2)
            solution2[j].Bin_Fullness = 0;
            solution2[j].L.free_linked_list();
        end

    end
    
    while(j < solution[number_items + 1].Bin_Fullness)
        solution[j].Bin_Fullness = 0;
        solution[j+=1].L.free_linked_list();
    end

    solution[number_items].Bin_Fullness = solution2[number_items].Bin_Fullness;
    solution[number_items + 1].Bin_Fullness = solution2[number_items + 1].Bin_Fullness;
    solution[number_items + 2].Bin_Fullness = solution2[number_items + 2].Bin_Fullness;
    solution[number_items + 3].Bin_Fullness = solution2[number_items + 3].Bin_Fullness;
    solution[number_items + 4].Bin_Fullness = solution2[number_items + 4].Bin_Fullness;
   
    if(delete_solution2)
        solution2[number_items].Bin_Fullness = 0;
        solution2[number_items + 1].Bin_Fullness = 0;
        solution2[number_items + 2].Bin_Fullness = 0;
        solution2[number_items + 3].Bin_Fullness = 0;
        solution2[number_items + 4].Bin_Fullness = 0;
    end

end



"""
/************************************************************************************************************************
 To free the memory of the individuals of the population                                            *
************************************************************************************************************************/
"""
function Clean_population()

  # long int  i,
  #     j;

    for(i = 0; i < P_size; i+=1)
        for(j = 0; j < number_items + 5; j+=1)    
            population[i][j].L.free_linked_list();
            population[i][j].Bin_Fullness = 0;
            children[i][j].L.free_linked_list();
            children[i][j].Bin_Fullness = 0;
        end
   end
end



"""
/************************************************************************************************************************
 To check if any of the items of the current bin is already in the solution                             *
 Input:                                                                                                                 *
  The position in the population of the solution: individual                                        *
  A new bin that could be added to the solution: bin                                              *
   An array that indicates the items that are already in the solution: items                                            *
 Output:                                                                                                                *
  (1) when none of the items in the current bin is already in the solution                                *
   (0) otherwise                                                                      *
************************************************************************************************************************/
"""
function Used_Items(long int individual, long int bin, long int items[])
  # long int   item,
  #     i,
  #     counter = 0;
  #  node *p;

    p = population[individual][bin].L.first;
    while(p != NULL)
        item = p->data;
        p = p->next;
        if(items [item] != 1)
            items_auxiliary[counter+=1] = item;
            items[item] = 1;
        else   
            for(i = 0; i < counter ; i+=1)
                items[items_auxiliary[i]] = 0;
            end
            return 0;
        end
    end
    
    return 1;
end



"""
/************************************************************************************************************************
 To put together all the used bins of the solution                                                *
 Input:                                                                                                                 *
  The position in the population of the solution: individual                                        *
************************************************************************************************************************/
"""
function Adjust_Solution(long int individual)

  # long int  i = 0,
  #     j = 0,
  #        k;
    while(population[individual][i].Bin_Fullness > 0)
        i+=1;
    end
    
    for(j = i, k = i; j < number_items; j+=1, k+=1)
        
        if(j < population[individual][number_items + 1].Bin_Fullness) 
            
            while(population[individual][k].Bin_Fullness == 0)
                k+=1;
            end
            
            population[individual][j].L.first = NULL;
            population[individual][j].L.last = NULL;
            population[individual][j].Bin_Fullness = population[individual][k].Bin_Fullness;
            population[individual][j].L.get_linked_list(population[individual][k].L);
        else
            population[individual][j].Bin_Fullness = 0;
            population[individual][j].L.first = NULL;
            population[individual][j].L.last = NULL;
            population[individual][j].L.num = 0;
        end
    end
end
