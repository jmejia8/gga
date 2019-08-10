
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

    file::String
    nameC::String
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
    file="",
    nameC=""
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
        seed_permutation,file,nameC)
    
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
    status.save_bestSolution = configuration[11]==1

    seed_permutation = status.seed;
    seed_emptybin = status.seed;
    
    status.ordered_population = collect(1:status.P_size);
    status.random_individuals = collect(1:status.P_size);
    status.best_individuals = collect(1:status.P_size);

    status.is_optimal_solution = false;
    status.generation = 0;
    status.repeated_fitness = 0

    seed!(status.seed)

    status
end
