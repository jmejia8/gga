"""
/************************************************************************************************************************
 To produce a small modification in a solution.                                                 *
 Input:                                                                                                         *
  The position in the population of the solution to mutate: individual                                *
  The rate of change to calculate the number of bins to eliminate: k                                  *
  A value that indicates if the solution was cloned: is_cloned                                      *
************************************************************************************************************************/
"""
function Adaptive_Mutation_RP( individual, float k, int is_cloned)

  # 
  #        number_bins,
  #        i,
  #        lightest_bin = 0,
  #        number_free_items = 0,
  #        free_items[ATTRIBUTES] = {0},
  #        ordered_BinFullness[ATTRIBUTES] = {0};
  #  node *p;

    for(i = 1:population[individual].number_of_bins)
        ordered_BinFullness[i] = i;
    end
   
    if(is_cloned)
        Sort_Random(ordered_BinFullness,0, population[individual].number_of_bins);
    end

    Sort_Ascending_BinFullness(ordered_BinFullness, individual);
    i = 1;
    
    while(population[individual][ordered_BinFullness[i]].Bin_Fullness < bin_capacity && i < population[individual].number_of_bins)
        i+=1;
    end    
    
    _p_ = 1 / (float)(k);
    number_bins = ceil(i*((2 - i/population[individual].number_of_bins) / pow(i,_p_))*(1 - ((double)get_rand(&seed_emptybin,()ceil((1/pow(i,_p_))*100))/100)));
    
    for i = 1:number_bins
    
        p = population[individual][ordered_BinFullness[lightest_bin]].w.first;

        while(p != NULL)
            free_items[number_free_items+=1] = p->data;
            p = p->next;
        end

        population[individual][ordered_BinFullness[lightest_bin]].w.free_linked_list();
        population[individual][ordered_BinFullness[lightest_bin]].Bin_Fullness = 0;
        lightest_bin+=1;
    end
 
    population[individual].number_of_bins -= number_bins;
    number_bins = population[individual].number_of_bins;
 
    Adjust_Solution(individual);
    RP(individual, number_bins, free_items, number_free_items);
 
    population[individual][number_items+1].Bin_Fullness = number_bins;
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
function RP( individual,  &b,  F[],  number_free_items)

  # 
  #     i,
  #        k,
  #        k2,
  #        ban,
  #        sum = 0,
  #        total_free = 0,
  #        ordered_BinFullness[ATTRIBUTES] = {0},
  #     *new_free_items = new  [2];

   # node   *p,
   #    *s,
   #       *aux;

   higher_weight = weight[F[0]];
   lighter_weight = weight[F[0]];
   bin_i = b;
   population[individual].fitness = 0;
   population[individual].fully_bins = 0;
   population[individual].highest_avaliable = bin_capacity;

    for (i = 0; i < b; i+=1 )
        ordered_BinFullness[i] = i;
    end
    
    Sort_Random(ordered_BinFullness,0, b);
    Sort_Random(F,0, number_free_items);

    for (i = 0; i < b; i+=1)
        sum = population[individual][ordered_BinFullness[i]].Bin_Fullness;

        p = population[individual][ordered_BinFullness[i]].w.first;
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
                            
                            if(population[individual][ordered_BinFullness[i]].w.last == s)
                                population[individual][ordered_BinFullness[i]].w.last = aux;
                            end
                            
                            population[individual][ordered_BinFullness[i]].w.num-=1;
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
                            
                            if(population[individual][ordered_BinFullness[i]].w.last == s)
                                population[individual][ordered_BinFullness[i]].w.last = aux;
                            end

                            population[individual][ordered_BinFullness[i]].w.num-=1;
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
      
        if(population[individual][ordered_BinFullness[i]].Bin_Fullness < population[individual].highest_avaliable)
            population[individual].highest_avaliable = population[individual][ordered_BinFullness[i]].Bin_Fullness;
        end
        
        if(population[individual][ordered_BinFullness[i]].Bin_Fullness == bin_capacity)
            population[individual].fully_bins+=1;
        elseif(population[individual][ordered_BinFullness[i]].Bin_Fullness + weight[ordered_weight[number_items-1]] <= bin_capacity)
            if(ordered_BinFullness[i] < bin_i)
                bin_i = ordered_BinFullness[i];
            end
        end
    end

    for(i = 0; i < bin_i; i+=1)
        population[individual].fitness += pow((population[individual].bins[i].Bin_Fullness / bin_capacity), 2);
    end
   
    free(new_free_items);
    number_free_items += total_free;

    if(higher_.weight < 5*bin_capacity)
        Sort_Random(F,0, number_free_items);
    else
        Sort_Descending_Weights(F, number_free_items);
        lighter_weight = weight[F[number_free_items-1]];
    end

    if(lighter_weight > bin_capacity - population[individual].highest_avaliable)
        for(i = bin_i; i < b; i+=1)
            population[individual].fitness += pow((population[individual].bins[i].Bin_Fullness / bin_capacity), 2);
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
 To find the solution with the highest fitness of the population and update the global_best_solution              *
************************************************************************************************************************/
"""
function Find_Best_Solution()

  #   i,
  #     best_individual = 0;
    for(i = 0; i < P_size; i+=1)
        if(population[i].fitness > population[best_individual].fitness)
            best_individual = i;
        end
    end

    if(generation + 1 > 1)
        if(population[best_individual].fitness > global_best_solution.fitness)
            global_best_solution = Copy_Solution( population[best_individual], 0);
        end
    else
        global_best_solution = Copy_Solution( population[best_individual], 0);
    end
end



"""
/************************************************************************************************************************
 To sort the individuals of the population in ascending order of their fitness                              *
************************************************************************************************************************/
"""
function Sort_Ascending_IndividualsFitness()

  #  i,
  #      k = P_size - 1,
  #      i2 = 0,
  #    aux,
  #      ban = 1;

    while(ban)
        ban = 0;
        for(i = i2; i < k; i+=1)
            if(population[ordered_population[i]].fitness > population[ordered_population[i+1]].fitness)
                aux = ordered_population[i];
                ordered_population[i] = ordered_population[i+1];
                ordered_population[i+1] = aux;
                ban = 1;
            else if(population[ordered_population[i]].fitness == population[ordered_population[i+1]].fitness)
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
function Sort_Ascending_BinFullness( ordered_BinFullness[],  individual)
  #   m,
  #     k,
  #        temporary_variable,
  #        ban = 1;

    k = population[individual].number_of_bins - 1;
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
 To copy solution2 in solution                                                            *
 Input:                                                                                                                 *
  A solution to save the copied solution: solution                                              *                                                                                                                 *
  The solution to be copied: solution2                                                        *
  A value that indicates if the copied solution must be deleted: delete_solution2                         *
************************************************************************************************************************/
"""
function Copy_Solution(solution2, status, delete_solution2)
    solution = deepcopy(solution2)
    solution.generation = status.generation;

    return solution

end



"""
/************************************************************************************************************************
 To free the memory of the individuals of the population                                            *
************************************************************************************************************************/
"""
# function Clean_population()

#     for(i = 0; i < P_size; i+=1)
#         for(j = 0; j < number_items + 5; j+=1)    
#             population[i][j].w.free_linked_list();
#             population[i][j].Bin_Fullness = 0;
#             children[i][j].w.free_linked_list();
#             children[i][j].Bin_Fullness = 0;
#         end
#    end
# end



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
function Used_Items( individual,  bin,  items[])
  #    item,
  #     i,
  #     counter = 0;
  #  node *p;

    p = population[individual][bin].w.first;
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
function Adjust_Solution( individual)

  #   i = 0,
  #     j = 0,
  #        k;
    while(population[individual].bins[i].Bin_Fullness > 0)
        i+=1;
    end
    
    for(j = i, k = i; j < number_items; j+=1, k+=1)
        
        if(j < population[individual].number_of_bins) 
            
            while(population[individual][k].Bin_Fullness == 0)
                k+=1;
            end
            
            population[individual][j].w.first = NULL;
            population[individual][j].w.last = NULL;
            population[individual][j].Bin_Fullness = population[individual][k].Bin_Fullness;
            population[individual][j].w.get_linked_list(population[individual][k].w);
        else
            population[individual][j].Bin_Fullness = 0;
            population[individual][j].w.first = NULL;
            population[individual][j].w.last = NULL;
            population[individual][j].w.num = 0;
        end
    end
end
