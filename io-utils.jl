"""
/************************************************************************************************************************
 To print the performance of the procedure on a BPP instance in a data file                             *
************************************************************************************************************************/
"""
function WriteOutput(status)
    output = fopen(status.nameC, "a");
    str = @sprintf("\n%s \t %d \t %d \t %d \t %f", file, L2, global_best_solution.fitness, generation, TotalTime)
    write(output, str);
    save_bestSolution && sendtofile(global_best_solution);
end

"""
/************************************************************************************************************************
 To print the global best solution in a data file                                               *
************************************************************************************************************************/
"""
function sendtofile(best)
    string1 = ""
    fil = ""
    aux = ""

    accumulated = 0;
    ban = 1,
    item = 0,
    position = 0,
    bins[ATTRIBUTES] = {0},
    n_bins = best[number_items + 1].Bin_Fullness;

    binError = -1,
    banError = 0;



    string(fil, "Details_GGA-CGT/GGA-CGT_S_(");
    string(string1, file);
    # itoa(conf, aux, 10);
    @sprintf(aux,"%ld", conf);
    string(fil,aux);
    string(fil,")_");
    string(fil,string1);
    try
        output = open(fil, "w+")
    catch
        return 
    end

    write(output,"Instance:\t%s\n", file);
    write(output,"Number of items:\t%ld\n", number_items);
    write(output,"Bin capacity:\t%ld\n", bin_capacity);
    write(output,"L2:\t%ld\n", L2);
    write(output,"\n****************************GGA-CGT global best solution******************************\n");
    write(output,"Number of bins:\n%ld\n", n_bins);
    write(output,"Fitness:\n%f\n", best[number_items].Bin_Fullness);
    write(output,"Optimal order of the weights:\n");

  # for (bin = 1:n_bins)
  # {  bins[bin] = 0;
  #   p = best[bin].L.first;
  #   while(true)
  #   {  if(p == NULL)
  #         break;
  #        item = p->data;
  #        p = p->next;
  #     bins[bin] += weight[item];
  #        accumulated += weight[item];
  #        fprintf(output, "%ld\n",weight[item]);
  #     if(bins[bin] > bin_capacity)
  #     {
  #       printf("ERROR the capacity of bin %ld was exceeded", bin);
  #           binError = bin;
  #       getchar();
  #           banError = 1;
  #     }

  #   }
  # }

  #  if(accumulated != total_accumulated_weight)
  #  {   printf("ERROR inconsistent sum of weights");
  #      getchar();

  #  }

  # fprintf(output,"\nDetailed solution:");
  # for (j=0; j < n_bins; j++)
  # {
  #   if()bins[j] > )bin_capacity)
  #       fprintf(output, " \n ********************ERROR the capacity of the bin was exceeded******************");
  #   fprintf(output, "\n\nBIN %ld\nFullness: %ld Gap: %ld\nStored items:\t ",j + 1, bins[j], bin_capacity-bins[j]);
  #     p = best[j].L.first;
  #   for (position=0; ; position++)
  #   {  if(p == NULL)
  #         break;
  #     item = p->data;
  #     p = p->next;
  #     fprintf(output, "[Item: %ld, Weight: %ld]\t", item + 1,weight[item]);
  #   }

  # }

  # fclose(output);

  #  if(banError)
  #   exit(1);
end
