function WriteOutput(status)
    open(status.nameC, "a") do output
      str = @sprintf("\n%s \t %d \t %d \t %d \t %f", status.file, status.L2, status.global_best_solution.fitness, status.generation, status.TotalTime)
      write(output, str)
    end
    status.save_bestSolution && sendtofile(status.global_best_solution, status);
end

function sendtofile(best, status)
    fil = string("Details_GGA-CGT/GGA-CGT_S_($(status.conf))_", status.file);

    open(fil, "w+") do output
      write(output,"Instance:\t$(status.file)\n");
      write(output,"Number of items:\t$(status.number_items)\n");
      write(output,"Bin capacity:\t$(status.bin_capacity)\n");
      write(output,"L2:\t$(status.L2)\n");
      write(output,"\n****************************GGA-CGT global best solution******************************\n");
      write(output,"Number of bins:\n$(best.number_of_bins)\n");
      write(output,"Fitness:\n$(best.fitness)\n");
      write(output,"Optimal order of the weights:\n");

      for bin = best.bins
        write(output, sprint(println, status.weight[bin.w]));
      end
    end

end
