#!/contrib/projects/julia/julia
using Distributions

mutable struct Lineage
  fitness::Float64
  size::Int64
  state::Vector{Int64}
  ploidy::Int64
end

function average_fitness(population::Dict{Array{Int64,1}, Lineage})
#calculate average fitness of the population. also outputs population size
  popN = Int64
  popN = 0
  popW = Float64
  popW = 0
  for line in values(population)
    popN+=line.size
    popW+=line.size * line.fitness
  end
  popW = popW/popN
  return popW, popN
end

function pop_size(population::Dict{Array{Int64,1}, Lineage})
  popN = Int64
  popN = 0
  for line in values(population)
    popN+=line.size
  end
  return popN
end


function hap_freq(population::Dict{Array{Int64,1}, Lineage})
#calculates haploid frequency
  popN = Int64
  popN = 0
  hapN = Int64
  hapN = 0
  for line in values(population)
    popN+=line.size
    if line.ploidy == 1 hapN+=line.size end
  end
  hapF = Float64
  hapF = hapN/popN
  return hapF
end

function ben_freq(population::Dict{Array{Int64,1}, Lineage})
#calculates beneficial frequency
  popN = Int64
  popN = 0
  benN = Int64
  benN = 0
  avgBen = Int64
  avgBen = 0
  for line in values(population)
    popN+=line.size
    if line.state[2]>0 benN+=line.size end
    avgBen+=line.state[2]*line.size
  end
  benF = Float64
  benF = benN/popN
  return benF, avgBen/popN, benN
end

function hap_ben(population::Dict{Array{Int64,1}, Lineage})
#checks if the haploid population carries at least one beneficial mutation

    check = false
    for line in values(population)
        if line.ploidy == 1
            if line.state[2]>0 check = true end
        end
    end
    return check
end

function estab_count(population::Dict{Array{Int64,1}, Lineage}, popw::Float64)
#checks for established beneficial lineages (frequency above 1/s)
  established_lineages = 0
  for line in values(population)
      if line.fitness > 1.0
          if line.size > 1/(line.fitness-1)
              established_lineages+=1
          end
      end
  end
  return established_lineages
end

function fitness_calc(state::Array{Int64,1}, sb::Float64,hb::Float64, sd::Float64,hd::Float64)
  w = Float64
  if state[1] == 2
  w=1.0+state[2]*(sb*hb)-state[3]*(sd*hd)
  else
  w=1.0+state[2]*(sb)-state[3]*(sd)
  end
  return w
end

function mutate_population(population::Dict{Array{Int64,1}, Lineage}, Ub::Float64, sb::Float64, Ud::Float64, sd::Float64, hb::Float64, hd::Float64)
  new_population = Dict{Array{Int64,1}, Lineage}()
  #initialize new population dictionary
  for (key, line) in population
    bmutations = rand(Poisson(line.ploidy * Ub * line.size))
    dmutations = rand(Poisson(line.ploidy * Ud * line.size))
      #for each new mutation, copies the state, adds mutation and generates the key.
      #check if it is already in the new_population. if it is, add 1 individual to it, if it's not - make another lineage with
      #the key and size 1
    for i = 1:bmutations
      #assign beneficial mutations
      newstate = copy(line.state)
      newstate[2]+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit::Float64 = fitness_calc(newstate,sb,hb, sd,hd)
        new_population[newstate] = Lineage(newfit, 1, newstate, line.ploidy)
      end
    end

    for i = 1:dmutations
      #assign deleterious mutations
      newstate = copy(line.state)
      newstate[3]+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit::Float64 = fitness_calc(newstate,sb,hb, sd,hd)
        new_population[newstate] = Lineage(newfit, 1, newstate, line.ploidy)
      end
    end
    new_size = line.size - bmutations - dmutations
    if  new_size > 0
      if key in keys(new_population) new_population[key].size+=new_size
      else new_population[key] = Lineage(line.fitness, new_size, line.state, line.ploidy)
      end
    end
  end
  return new_population
end

function wright_fisher_reproduction(population::Dict{Array{Int64,1}, Lineage}, N0::Int64, Nnew::Int64, popw::Float64)
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[] #initialize lineage liste
  for (key,line) in population
    push!(proby_list, line.size/N0 * line.fitness/popw)
    #representation of a lineage in the next generations depends on its frequency and relative fitness
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(Nnew, proby_list))
  new_population = Dict{Array{Int64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].ploidy)
    end
  end
  return new_population
end

function simulate(pop_Ni::Int64, hb::Float64, hd::Float64)
  #parameters are defined here
  init_hap_N = 1
  sb = 0.1
  sd = 0.1
  Ub = 0.0001
  Ud = 0.01
  population = Dict{Array{Int64,1}, Lineage}()
  hap_state = [1,0,0]  #lineage state. 1: ploidy, 2: ben count, 3: del count
  population[hap_state] = Lineage(1.0,init_hap_N,hap_state,1)
  dip_state = [2,0,0]
  population[dip_state] = Lineage(1.0,pop_Ni-init_hap_N,dip_state,2)
  hap_f = hap_freq(population)
  popw, popN = average_fitness(population)
  generations = Int64
  generations = 0 #initialize generations counter
  bens_in_haps=false
  while 0.0<hap_f<1.0
    popw, pop_N = average_fitness(population)
    population = wright_fisher_reproduction(population, pop_N, pop_Ni, popw)
    population = mutate_population(population, Ub, sb, Ud, sd, hb, hd)
    if bens_in_haps == false bens_in_haps = hap_ben(population) end
    hap_f = hap_freq(population)
    popw, pop_N = average_fitness(population)
    generations+=1
  end
  ben_state_ret = ben_freq(population)
  if ben_state_ret[3] > 1/(sb*hb) est=1
  else est=0
  end
  return hap_f ==1.0, generations, popw, ben_state_ret[1],ben_state_ret[2],est, bens_in_haps
end

job_id = ARGS[1]
hd=0.0
for hb in (0.0, 0.25, 0.5, 0.75, 1.0)
    for n in [10, 14, 20, 29, 42, 61, 88, 127, 183, 263, 379, 545, 784, 1128, 1623, 2335, 3359, 4832, 6951, 10000]
        run_time= Int64[] #run_time list records all times of haploid fixation
        est_count= Int64[]
        fitness_endW = Float64[] #just haploid hitchhikings
        ben_hap = Int64[]
        for run = 1:1000
            output = simulate(n, hb, hd)
            if output[1]
              push!(run_time, output[2])
              push!(fitness_endW, output[3])
              if output[7] push!(ben_hap, 2)
              else push!(ben_hap, 1)
              end
            else
              push!(run_time, 0)
              push!(est_count, output[6])
              if output[7] push!(ben_hap, 3)
              else push!(ben_hap, 0)
              end
            end
        end
        outfile = open(string(job_id,"run_time.csv"), "a")
        write(outfile, join(run_time, ","), "\n")
        close(outfile)
        outfile = open(string(job_id,"fitness_endW.csv"), "a")
        write(outfile, join(fitness_endW, ","), "\n")
        close(outfile)
        outfile = open(string(job_id,"haps_bens.csv"), "a")
        write(outfile, join(ben_hap, ","), "\n")
        close(outfile)
        outfile = open(string(job_id,"est_count.csv"), "a")
        write(outfile, join(est_count, ","), "\n")
        close(outfile)
    end
end
