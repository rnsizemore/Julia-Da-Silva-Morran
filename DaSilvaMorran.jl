import Plots
using Statistics
# parameters
IL = 30             # num. of host interaction loci
L = IL + 1          # num. of host loci
initfr = 0.2        # inital freq. of recomb. allele
standing = true     # create standing genetic variation
s = 0.8             # max selection coeff.
u = 10^-4          # mutation rate per locus
N = [1500, 750, 750]            # host pop.
r = [0.5, 0.5, 0.1]             # recombination rate
ngens = 30          # num. of generations
npgens = 1          # parasite gen. per host gen.
nreps = 30           # num. of simulations to run
nscens = 3

# variables
hosts = zeros(Int, N[1], L)                                    # holds value for each locus in each individual
hostsnew = zeros(Int, N[1], L)                                 # holds hosts value for next generation
temp = Array{Int}(undef, L)                                 # temp holder for an individuals genes
host0 = zeros(L)                                            # freq. of mutant allele in population at each locus
host1 = zeros(L)                                            # freq. of wild type allele in population at each locus
para0 = zeros(L-1)                                            # freq. of parasite mutant allele
para1 = zeros(L-1)                                            # freq. of parasite wild type allele
w0 = Array{Float64}(undef, L-1)                                 # host allele fitness for mutant alleles
w1 = Array{Float64}(undef, L-1)                                 # host allele fitness for wild type alleles
wp0 = Array{Float64}(undef, L-1)                                # parasite mutant allele fitnesses
wp1 = Array{Float64}(undef, L-1)                                # parasite wild type allele fitnesses
wp_bar = Array{Float64}(undef, L-1)                             # parasite mean fitness
w = Array{Float64}(undef, N[1])                                    # individual host fitnesses
fix = false                                                 # flag for fixation of recombination allele
nfix = 0                                                    # number of simulations where recombination allele fixates
gfix = zeros(Int, nreps)                                    # stores the generation where recombination allele comes to fixation
frecomb = Array{Float64}(undef, nreps)                          # freq. of recombination allele at end of simulation run
genrecomb = Array{Float64}(undef, ngens, nscens)                        # average freq. of recombination allele at each generation over the replicates
genrecombarr = Array{Float64}(undef, ngens, nreps, nscens)              # stores the array of values of recomb. freq. at each generation
genrecombstd = Array{Float64}(undef, ngens, nscens)                     # stores standard deviation of replicates recomb. freq.
hostMutantAlleleFreq = Array{Float64}(undef, nreps, ngens, L)   # host mutant allele freq. for each generation of a simulation run
meanf = Array{Float64}(undef, ngens, L)                         # mean host mutant allele freq. over all simulation runs
iFixCount = zeros(Int, nreps, L-1)                          # number of simulation runs where each locus has fixation
iPrevFix = Array{Int}(undef, L-1)                           # the last fixed allele for each host interaction locus
rep = 1             # current simulation count
gen = 1             # current generation count 
scen = 1            # current scenario
pgen = 1            # current parasite generation count
printc = false      # flag to print information
random_genomes = rand(Float64, Int(N[1]/2))    # random number for inital freq of mutant allele

for scen in 1:3
    for rep = 1:nreps
        global hostsnew = zeros(Int, N[scen], L)
        # logic for printing when running many repetitions of the simulation
        if(nreps <= 10) # print every run if doing less than 10
            global printc = true
        else
            if(nreps <= 100)
                if(mod(rep, 10) == 0)   # print every 10 runs if doing less than 100
                    global printc = true
                end
            else
                if(mod(rep, 100) == 0)  # print every 100 runs if doing more
                    global printc = true
                end
            end
        end
        if(printc)
            print("rep ", rep, "\n")
        end

        global hosts = zeros(Int, N[scen], L)
        for i = 1:Int(N[scen]*initfr) # put inital recomb allele in with initfr freq.
            hosts[i, 1] = 1
        end

        if(standing) # create standing variation
            for i = 2:L
                for j = 1:Int(N[scen]/2)
                    hosts[Int(round(random_genomes[j]*N[scen], RoundDown)) + 1, i] = 1  # make random_genomes proportion of individuals have mutant alleles
                end
            end
        end

        for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
            host1[i] = sum(hosts[:, i])/N[scen]
            host0[i] = 1.0 - host1[i]
        end

        # intitalize wild type allele at fixation in parasite population
        para0 .= 1
        para1 .= 0

        global fix = false                 # reset fixation of recomb allele
        global iPrevFix = zeros(Int, L-1)  # intitalize last fixed allele as 0

        for gen = 1:ngens # loop through generations
            #print("\n Gen #", gen, "\n")
            #print("Start of gen: ", sum(hosts[:, 1])/N[scen], "\n")
            for i = 1:L-1
                w0[i] = 1 - s*para0[i]  # wild type allele fitness
                w1[i] = 1 - s*para1[i]  # mutant allele fitness
            end
            global w = ones(N[scen])          # intitalize individual fitnesses
            for i = 1:N[scen]
                for j = 1:L-1
                    w[i] = w[i] * ( w0[j] * ( 1 - hosts[i, j+1] ) + w1[j] * hosts[i, j+1] ) # calculate individual fitnesses
                end
            end

            for i = 1:L-1   # calculate parasite fitnesses
                wp0[i] = 1.0 - s*host1[i+1]
                wp1[i] = 1.0 - s*host0[i+1]
            end

            # host reproduction and selection
            wMax = maximum(w)
            for i = 1:N[scen]
                w[i] = w[i] / wMax  # calculate relative fitnesses
            end
            #print("\n Wild type fitness: ", w0[1], "\n")
            #print(" Wild type parasite fitness: ", wp0[1], "\n \n")
            for i = 1:N[scen]
                while(true)
                    j = rand(1:N[scen])
                    if(rand() <= w[j])
                        hostsnew[i, :] = hosts[j, :]
                        break
                    end
                end
            end
            
            global hosts = hostsnew
            #print("After reproduction: ", sum(hosts[:, 1])/N[scen], "\n")
            # mutation
            for i = 1:N[scen]
                for j = 1:L # loop through loci in all individuals
                    if(rand() <= u) # if random num is less than mutation rate
                        hosts[i, j] = abs(hosts[i, j] - 1)  # flip the allele of the locus
                    end
                end
            end
            #print("After mutation: ", sum(hosts[:, 1])/N[scen], "\n")
            # recombination
            for i = 1:2:N[scen]
                if(hosts[i, 1] == 1 || hosts[i + 1, 1] == 1) # loop through pairs of hosts, checking for recombination allele
                    for j = 1:L-1
                        if(rand() <= r[scen]) # 0.5 chance at each locus to swap between pairs
                            global temp = hosts[i, :]
                            hosts[i, j+1:L] = hosts[i+1, j+1:L]
                            hosts[i + 1, j+1:L] = temp[j+1:L]
                        end
                    end
                end
            end
            #print("After recombination: ", sum(hosts[:, 1])/N, "\n")
            # parasites
            for i = 1:L-1   # calculate parasite fitness
                wp0[i] = 1 - s*host1[i+1]
                wp1[i] = 1 - s*host0[i+1]
            end
            for pgen = 1:npgens
                for i = 1:L-1
                    wp_bar[i] = wp0[i] * para0[i] + wp1[i] * para1[i]  # mean parasite fitness
                    para0[i] = para0[i] * wp0[i] / wp_bar[i]         # freq. of parasite wild type alleles
                    para1[i] = 1.0 - para0[i]                     # freq. of parasite mutant alleles

                    para0[i] = para0[i] * (1.0 - u) + para1[i] * u # wild type allele freq. after mutation
                    para1[i] = 1.0 - para0[i]                     # mutant type allele freq. after mutation
                end
            end

            # fixation
            for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
                host1[i] = sum(hosts[:, i])/N[scen]
                host0[i] = 1.0 - host1[i]
            end
            hostMutantAlleleFreq[rep, gen, :] .= host1  # store mutant allele freq.
            frecomb[rep] = host1[1]     # store ending freq. of recombination allele
            
            if(!fix && host1[1] > 0.99)
                global fix = true
                global nfix += 1
                gfix[rep] = gen
                #break
            end
            
            for i = 1:L-1
                if(host1[i+1] > 0.99 && iPrevFix == 0)
                    iFixCount[rep, :] = iFixCount[rep, :] .+ 1
                    iPrevFix = 1
                end
                if(host0[i+1] > 0.99 && iPrevFix == 1)
                    iFixCount[rep, :] = iFixCount[rep, :] .+ 1
                    iPrevFix = 0
                end
            end
            genrecomb[gen, scen] += host1[1]/nreps         # add this generation's recomb allele freq. to genrecomb
            genrecombarr[gen, rep, scen] = host1[1]        # add this generation's recomb allele freq. to genrecombarr
            #print("End of gen: ", sum(hosts[:, 1])/N, "\n")
        end
    end
    for i = 1:ngens
        genrecombstd[i] = std(genrecombarr[i, :, scen])
    end
end
println("\nFixations at Host Interaction Loci")
println("No. of fixations/locus/gen: ",  sum(iFixCount) / (IL * nreps) / ngens)
println("Apparent no. of fixations/locus/gen: ", count(hostMutantAlleleFreq[:, ngens, 2:L] .> 0.99) / (nreps * IL) / ngens, "\n")

println("Fixation of Recombination Allele")
println("Prop. of replicates with fixation: ", nfix / nreps)
if(nfix > 0)
    println("Mean no. of generations to fixation: ", sum(gfix) / nfix)
end
println("Mean freq of recomb allele: ", sum(frecomb) / nreps)

meanf = sum(hostMutantAlleleFreq[:])
#@show genrecombarr[ngens, :]
#@show genrecombstd

using Plots
testPlot = plot(1:ngens, genrecomb[:, 1])
plot!(genrecomb[:, 2])
plot!(genrecomb[:, 3])
ylims!(0, 1)
display(testPlot)
readline()