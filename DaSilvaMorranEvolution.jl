# evolution experiment, stable parasite present 

using Plots
using Statistics
using CSV
using DataFrames
# parameters
IL = 30             # num. of host interaction loci
L = IL + 1          # num. of host loci
initfr = 0.2        # inital freq. of recomb. allele
standing = true     # create standing genetic variation
s = 0.8             # max selection coeff.
u = 10^-4          # mutation rate per locus
NArray = [1500, 750, 750]            # host pop.
rArray = [0.5, 0.5, 0.1]             # recombination rate
ngens = 30          # num. of generations
npgens = 1          # parasite gen. per host gen.
nreps = 100           # num. of simulations to run
nscens = 3

# variables
hosts = zeros(Int, NArray[1], L)                                    # holds value for each locus in each individual
hostsnew = zeros(Int, NArray[1], L)                                 # holds hosts value for next generation
temp = Array{Int}(undef, L)                                 # temp holder for an individuals genes
host0 = zeros(L)                                            # freq. of mutant allele in population at each locus
host1 = zeros(L)                                            # freq. of wild type allele in population at each locus
para0 = zeros(L-1)                                            # freq. of parasite mutant allele
para1 = zeros(L-1)                                            # freq. of parasite wild type allele
w0 = Array{Float64}(undef, L)                                 # host allele fitness for mutant alleles
w1 = Array{Float64}(undef, L)                                 # host allele fitness for wild type alleles
wp0 = Array{Float64}(undef, L-1)                                # parasite mutant allele fitnesses
wp1 = Array{Float64}(undef, L-1)                                # parasite wild type allele fitnesses
wp_bar = Array{Float64}(undef, L-1)                             # parasite mean fitness
w = Array{Float64}(undef, NArray[1])                                    # individual host fitnesses
fix = false                                                 # flag for fixation of recombination allele
nfix = 0                                                    # number of simulations where recombination allele fixates
gfix = zeros(Int, nreps)                                    # stores the generation where recombination allele comes to fixation
frecomb = Array{Float64}(undef, nreps)                          # freq. of recombination allele at end of simulation run
genrecomb = zeros(Float64, ngens+1, nscens)                        # average freq. of recombination allele at each generation over the replicates
genrecombarr = zeros(Float64, ngens+1, nreps, nscens)              # stores the array of values of recomb. freq. at each generation
genrecombstd = zeros(Float64, ngens+1, nscens)                     # stores standard deviation of replicates recomb. freq.
hostMutantAlleleFreq = Array{Float64}(undef, nreps, ngens+1, L)   # host mutant allele freq. for each generation of a simulation run
meanf = Array{Float64}(undef, ngens+1, L)                         # mean host mutant allele freq. over all simulation runs
iFixCount = zeros(Int, nreps, L-1)                          # number of simulation runs where each locus has fixation
iPrevFix = Array{Int}(undef, L)                           # the last fixed allele for each host interaction locus
rep = 1             # current simulation count
gen = 1             # current generation count 
scen = 1            # current scenario
pgen = 1            # current parasite generation count
printc = false      # flag to print information
#random_genomes = rand(Float64, Int(N[1]/2))    # random number for inital freq of mutant allele
#global 
printc = nreps <= 10 # print every run if doing less than 10
for scen in 1:3
    # local initialization steps for each scenario
    N = NArray[scen]
    r = rArray[scen]
    print("Scenario $scen\n")
    print("N=$N, r=$r\n")
    nfix = 0
    gfix = zeros(Int, nreps)
    hostsnew = zeros(Int, N, L)
    for rep = 1:nreps
        # logic for printing when running many repetitions of the simulation
        printc = true
        if nreps > 10 && nreps <= 100
            printc = mod(rep, 10) == 0 # print every 10 runs if doing less than 100
        elseif nreps > 100
            printc = mod(rep, 50) == 0 # print every 100 runs if doing more
        end
        if printc
            print("rep ", rep, "\n")
        end
        hosts = zeros(Int, N, L)
        for i = 1:Int(N*initfr) # put inital recomb allele in with initfr freq.
            hosts[i, 1] = 1
        end
        random_genomes = rand(Float64, Int(N/2))    # random number for inital freq of mutant allele
        if(standing) # create standing variation
            for i = 2:L
                random_genomes = rand(Float64, Int(N/2))
                for j = 1:Int(N/2)
                    hosts[Int(round(random_genomes[j]*N, RoundDown)) + 1, i] = 1  # make random_genomes proportion of individuals have mutant alleles
                end
            end
        end

        for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
            host1[i] = sum(hosts[:, i])/N
            host0[i] = 1.0 - host1[i]
        end
        # initialization of data collection structures
        gen = 1
        genrecomb[gen, scen] += host1[1]/nreps         # add generation 0's recomb allele freq. to genrecomb
        if genrecomb[gen,scen]>1 || genrecomb[gen,scen]<0
            error("Out of range value in genrecomb[$gen,$scen]")
        end
        genrecombarr[gen, rep, scen] = host1[1]        # add generation 0's recomb allele freq. to genrecombarr
        hostMutantAlleleFreq[rep, gen, :] = host1  # store mutant allele freq. at generation 0
        # intitalize wild type allele at fixation in parasite population
        para0 .= 0
        para1 .= 1

        fix = false                 # reset fixation of recomb allele
        iPrevFix = zeros(Int, L)  # intitalize last fixed allele as 0

        for gen = 1:ngens # loop through generations
            #print("\n Gen #", gen, "\n")
            #print("Start of gen: ", sum(hosts[:, 1])/N[scen], "\n")
            #for i = 2:L
                @. w0[2:L] = 1 - s*para0  # wild type allele fitness
                @. w1[2:L] = 1 - s*para1  # mutant allele fitness
            #end
            w = ones(N)          # intitalize individual fitnesses
            for i = 1:N
                for j = 2:L
                    w[i] = w[i] * ( w0[j] * ( 1 - hosts[i, j] ) + w1[j] * hosts[i, j] ) # calculate individual fitnesses
                end
            end

            #for i = 2:L   # calculate parasite fitnesses
                @. wp0 = 1.0 - s*host1[2:L]
                @. wp1 = 1.0 - s*host0[2:L]
            #end

            # host reproduction and selection
            wMax = maximum(w)
            for i = 1:N
                w[i] = w[i] / wMax  # calculate relative fitnesses
            end
            #print("\n Wild type fitness: ", w0[1], "\n")
            #print(" Wild type parasite fitness: ", wp0[1], "\n \n")
            for i = 1:N
                while(true)
                    j = rand(1:N)
                    if(rand() <= w[j])
                        hostsnew[i, :] = hosts[j, :]
                        break
                    end
                end
            end
            
            hosts = hostsnew
            #print("After reproduction: ", sum(hosts[:, 1])/N[scen], "\n")
            # mutation
            for i = 1:N
                for j = 1:L # loop through loci in all individuals
                    if(rand() <= u) # if random num is less than mutation rate
                        hosts[i, j] = abs(hosts[i, j] - 1)  # flip the allele of the locus
                    end
                end
            end
            #print("After mutation: ", sum(hosts[:, 1])/N[scen], "\n")
            # recombination
            for i = 1:2:N
                if(hosts[i, 1] == 1 || hosts[i + 1, 1] == 1) # loop through pairs of hosts, checking for recombination allele
                    for j = 1:L-1
                        if(rand() <= r) # 0.5 chance at each locus to swap between pairs
                            temp = hosts[i, :]
                            hosts[i, j+1:L] = hosts[i+1, j+1:L]
                            hosts[i + 1, j+1:L] = temp[j+1:L]
                        end
                    end
                end
            end
            #print("After recombination: ", sum(hosts[:, 1])/N, "\n")
            # parasites
            #for i = 1:L-1   # calculate parasite fitness
                @. wp0 = 1 - s*host1[2:L]
                @. wp1 = 1 - s*host0[2:L]
            #end
            #for pgen = 1:npgens
                #for i = 1:L-1
            #        @. wp_bar = wp0 * para0 + wp1 * para1  # mean parasite fitness
            #        @. para0 = para0 * wp0 / wp_bar         # freq. of parasite wild type alleles
            #        @. para1 = 1.0 - para0                     # freq. of parasite mutant alleles

            #        @. para0 = para0 * (1.0 - u) + para1 * u # wild type allele freq. after mutation
            #        @. para1 = 1.0 - para0                     # mutant type allele freq. after mutation
                #end
            #end

            # fixation
            for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
                host1[i] = sum(hosts[:, i])/N
                host0[i] = 1.0 - host1[i]
            end
            hostMutantAlleleFreq[rep, gen+1, :] = host1  # store mutant allele freq.
            frecomb[rep] = host1[1]     # store ending freq. of recombination allele
            
            if(!fix && host1[1] > 0.99)
                fix = true
                nfix += 1
                gfix[rep] = gen
                #break
            end
            
            for i = 2:L
                if(host1[i] > 0.99 && iPrevFix[i] == 0)
                    iFixCount[rep, :] = iFixCount[rep, :] .+ 1
                    iPrevFix[i] = 1
                end
                if(host0[i] > 0.99 && iPrevFix[i] == 1)
                    iFixCount[rep, :] = iFixCount[rep, :] .+ 1
                    iPrevFix[i] = 0
                end
            end
            genrecomb[gen+1, scen] += host1[1]/nreps         # add this generation's recomb allele freq. to genrecomb
            if genrecomb[gen+1,scen]>1 || genrecomb[gen+1,scen]<0
                error("Out of range value in genrecomb[$(gen+1),$scen]")
            end
            genrecombarr[gen+1, rep, scen] = host1[1]        # add this generation's recomb allele freq. to genrecombarr
            #print("End of gen: ", sum(hosts[:, 1])/N, "\n")
        end
    end
    for i = 0:ngens
        genrecombstd[i+1] = std(genrecombarr[i+1, :, scen])
    end
    println("\nFixations at Host Interaction Loci")
    println("No. of fixations/locus/gen: ",  sum(iFixCount) / (IL * nreps) / (ngens+1))
    println("Apparent no. of fixations/locus/gen: ", count(hostMutantAlleleFreq[:, ngens+1, 2:L] .> 0.99) / (nreps * IL) / (ngens+1), "\n")

    println("Fixation of Recombination Allele")
    println("Prop. of replicates with fixation: ", nfix / nreps)
    if(nfix > 0)
        println("Mean no. of generations to fixation: ", sum(gfix) / nfix)
    end
    println("Mean freq of recomb allele: ", sum(frecomb) / nreps)
    println("\n")
end
meanf = sum(hostMutantAlleleFreq[:])
#@show genrecombarr[ngens, :]
#@show genrecombstd
size(genrecomb)
labels = reshape(map((x,y)->"N=$x, r=$y",NArray,rArray),1,nscens)
gr()
testPlot = plot(0:ngens, genrecomb,label=labels,plot_title="$nreps reps per scenario",
                xlabel="Generation",ylabel="Frequency",ylims=(0,1));
#plot!(genrecomb[:, 2])
#plot!(genrecomb[:, 3])
#ylims!(testPlot,0, 1)
#xlabel!(testPlot,"Generation")
#ylabel!(testPlot,"Frequency")

display(testPlot)
#readline()
dir_name = "plots"
if !isdir(dir_name)
    mkpath(dir_name)
end
plot_file = "nreps$(nreps)_ngens$(ngens)_evo"
savefig(joinpath(dir_name, "$(plot_file).pdf"))

##################################################
# Now let's plot the original Fortran simulations 
#labels = reshape(map((x,y)->"N=$x, r=$y",NArray,rArray),1,nscens)
#orig_folder = "../doi_10_5061_dryad_6343j__v20170310/"
#file_names = map((x,y)->joinpath(orig_folder,"meanHostMutantAlleleFreq_N$(x)_r$(y).out"),NArray,rArray)
#generecombFort = Array{Float64}(undef,ngens+1,nscens)
#for ni in 1:nscens
#   generecombFort[2:ngens+1,ni] = CSV.read(file_names[ni],DataFrame,header=1,delim=' ',ignorerepeated=true).L001
#end
# let's add the initial generation (gen 0)
#generecombFort[1,:] = repeat([initfr],inner=nscens)
#size(generecombFort)
#gr()
#testPlot = plot(0:ngens, generecombFort,label=labels,plot_title="$nreps reps per scenario",
#                xlabel="Generation",ylabel="Frequency",ylims=(0,1));

#display(testPlot)
#dir_name = "plots"
#if !isdir(dir_name)
#    mkpath(dir_name)
#end
#plot_file = "orig_nreps$(nreps)_ngens$ngens"
#savefig(joinpath(dir_name, "$(plot_file).pdf"))



#readline()