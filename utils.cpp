
"""
/************************************************************************************************************************
 To read the data defining a BPP instance                                                     *
************************************************************************************************************************/
"""
void LoadData()
{
    char  string[300];
  long  k;
   long int ban = 0;
   long double bin_capacity1;
   long double total_accumulated_aux = 0;

    FILE  *data_file;

   string[0] = '\0';
    strcpy(string, file);
  if((data_file = fopen(string, "rt")) == NULL)
    {
      printf("\nThere is no data file ==> [%s]%c", string, 7);
      return 0;
    }
   printf("\nThe file is %s", string);
    fgets(string, 300, data_file);
    fscanf(data_file, "%ld\n", &number_items);
    bin_capacity = 0;
    fgets(string, 300, data_file);
    fscanf(data_file, "%Lf\n", &bin_capacity1);
    best_solution = 0;
    fgets(string, 300, data_file);
    fscanf(data_file, "%ld\n", &best_solution);
    fgets(string, 300, data_file);
    total_accumulated_weight = 0;
    for(k = 0; k < number_items; k++)
    {   fscanf(data_file, "%Lf", &weight1[k]);
      weight[k] = (long int)weight1[k];
      total_accumulated_weight = (total_accumulated_weight + weight[k]);
      total_accumulated_aux += weight1[k];
      if(ban == 0)
    {
        if(weight1[k] / weight[k] > 1)
          ban = 1;
      }
    }
   if(ban)
   {  total_accumulated_weight = 0;
    for(k = 0; k < number_items; k++)
    { weight[k] = (long int)(ceil(weight1[k]*bin_capacity1 - .5) );
         total_accumulated_weight = (total_accumulated_weight + weight[k]);
      }
      bin_capacity1 *= bin_capacity1;
   }
   bin_capacity = (long int)bin_capacity1;
   fclose(data_file);
   if(ban)
   {  if((long int)total_accumulated_weight != (long int)(ceil(total_accumulated_aux*sqrt(bin_capacity) - .5) ))
    { printf("\t Error loading weights");
      # getch();
        exit(1);
    }
   }
   return 1;
}



"""
/************************************************************************************************************************
 To print the performance of the procedure on a BPP instance in a data file                             *
************************************************************************************************************************/
"""
function WriteOutput()
{
   output = fopen(nameC, "a");
   fprintf(output, "\n%s \t %d \t %d \t %d \t %f", file, (int)L2, (int)global_best_solution[number_items + 1].Bin_Fullness, generation, TotalTime);
  if(save_bestSolution == 1)
    sendtofile(global_best_solution);
   fclose(output);
}



"""
/************************************************************************************************************************
 To print the global best solution in a data file                                               *
************************************************************************************************************************/
"""
function sendtofile(SOLUTION best[])
{
  char  string1[300],
      fil[300],
         aux[300];

  long double accumulated = 0;
  long int   bin,
    ban = 1,
    item = 0,
    position = 0,
      bins[ATTRIBUTES] = {0},
    n_bins = best[number_items + 1].Bin_Fullness;

   int  binError = -1,
         banError = 0;
   long int j;



    FILE *output;
   node *p;
    strcpy(fil, "Details_GGA-CGT/GGA-CGT_S_(");
    strcpy(string1, file);
   # itoa(conf, aux, 10);
   sprintf(aux,"%ld", conf);
   strcat(fil,aux);
   strcat(fil,")_");
  strcat(fil,string1);
  if((output = fopen(fil, "w+")) == NULL)
    {
      # printf("\nThere is no data file ==> [%s]%c", file, 7);
    getchar();
    exit(1);
    }
  fprintf(output,"Instance:\t%s\n", file);
  fprintf(output,"Number of items:\t%ld\n", number_items);
  fprintf(output,"Bin capacity:\t%ld\n", bin_capacity);
  fprintf(output,"L2:\t%ld\n", L2);
   fprintf(output,"\n****************************GGA-CGT global best solution******************************\n");
  fprintf(output,"Number of bins:\n%ld\n", n_bins);
    fprintf(output,"Fitness:\n%f\n", best[number_items].Bin_Fullness);
  fprintf(output,"Optimal order of the weights:\n");

  for (bin = 0; bin < n_bins; bin++)
  {  bins[bin] = 0;
    p = best[bin].L.first;
    while(true)
    {  if(p == NULL)
          break;
         item = p->data;
         p = p->next;
      bins[bin] += weight[item];
         accumulated += weight[item];
         fprintf(output, "%ld\n",weight[item]);
      if(bins[bin] > bin_capacity)
      {
        printf("ERROR the capacity of bin %ld was exceeded", bin);
            binError = bin;
        getchar();
            banError = 1;
      }

    }
  }

   if(accumulated != total_accumulated_weight)
   {   printf("ERROR inconsistent sum of weights");
       getchar();

   }

  fprintf(output,"\nDetailed solution:");
  for (j=0; j < n_bins; j++)
  {
    if((long int)bins[j] > (long int)bin_capacity)
        fprintf(output, " \n ********************ERROR the capacity of the bin was exceeded******************");
    fprintf(output, "\n\nBIN %ld\nFullness: %ld Gap: %ld\nStored items:\t ",j + 1, bins[j], bin_capacity-bins[j]);
      p = best[j].L.first;
    for (position=0; ; position++)
    {  if(p == NULL)
          break;
      item = p->data;
      p = p->next;
      fprintf(output, "[Item: %ld, Weight: %ld]\t", item + 1,weight[item]);
    }

  }

  fclose(output);

   if(banError)
    exit(1);
}
