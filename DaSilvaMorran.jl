import Plots
using Statistics
# parameters
IL = 30             # num. of host interaction loci
L = IL + 1          # num. of host loci
initfr = 0.2        # inital freq. of recomb. allele
standing = true     # create standing genetic variation
s = 0.8             # max selection coeff.
u = 0.0001          # mutation rate per locus
N = 750            # host pop.
r = 0.5             # recombination rate
ngens = 30          # num. of generations
npgens = 1          # parasite gen. per host gen.
nreps = 50           # num. of simulations to run

# variables
hosts = zeros(Int, N, L)                                    # holds value for each locus in each individual
hostsnew = zeros(Int, N, L)                                 # holds hosts value for next generation
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
w = Array{Float64}(undef, N)                                    # individual host fitnesses
fix = false                                                 # flag for fixation of recombination allele
nfix = 0                                                    # number of simulations where recombination allele fixates
gfix = zeros(Int, nreps)                                    # stores the generation where recombination allele comes to fixation
frecomb = Array{Float64}(undef, nreps)                          # freq. of recombination allele at end of simulation run
genrecomb = Array{Float64}(undef, ngens)                        # average freq. of recombination allele at each generation over the replicates
genrecombarr = Array{Float64}(undef, ngens, nreps)              # stores the array of values of recomb. freq. at each generation
hostMutantAlleleFreq = Array{Float64}(undef, nreps, ngens, L)   # host mutant allele freq. for each generation of a simulation run
meanf = Array{Float64}(undef, ngens, L)                         # mean host mutant allele freq. over all simulation runs
iFixCount = zeros(Int, nreps, L-1)                          # number of simulation runs where each locus has fixation
iPrevFix = Array{Int}(undef, L-1)                           # the last fixed allele for each host interaction locus
rep = 1             # current simulation count
gen = 1             # current generation count
pgen = 1            # current parasite generation count
printc = false      # flag to print information
random_genomes = rand(Float64, Int(N/2))    # random number for inital freq of mutant allele

for rep = 1:nreps
    # reinitalize things
    global hosts = zeros(Int, N, L)                                    # holds value for each locus in each individual
    global hostsnew = zeros(Int, N, L)                                 # holds hosts value for next generation
    global temp = Array{Int}(undef, L)                                 # temp holder for an individuals genes
    global host0 = zeros(L)                                            # freq. of mutant allele in population at each locus
    global host1 = zeros(L)                                            # freq. of wild type allele in population at each locus
    global para0 = zeros(L-1)                                            # freq. of parasite mutant allele
    global para1 = zeros(L-1)                                            # freq. of parasite wild type allele
    
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
    for i = 1:Int(N*initfr) # put inital recomb allele in with initfr freq.
        hosts[i, 1] = 1
    end

    if(standing) # create standing variation
        for i = 2:L
            for j = 1:Int(N/2)
                hosts[Int(round(random_genomes[j]*N, RoundUp)), i] = 1  # make random_genomes proportion of individuals have mutant alleles
            end
        end
    end

    for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
        tempSum = 0
        for j = 1:N
            tempSum += hosts[j, i]
        end
        host1[i] = tempSum/N
        host0[i] = 1 - host1[i]
    end

    # intitalize wild type allele at fixation in parasite population
    para0 .= 1
    para1 .= 0

    global fix = false                 # reset fixation of recomb allele
    global iPrevFix = zeros(Int, L-1)  # intitalize last fixed allele as 0

    for gen = 1:ngens # loop through generations
        for i = 1:L-1
            w0[i] = 1 - s*para0[i]  # wild type allele fitness
            w1[i] = 1 - s*para1[i]  # mutant allele fitness
        end
        global w = ones(N)          # intitalize individual fitnesses
        for i = 1:N
            for j = 1:L-1
                w[i] = w[i] * ( w0[j] * ( 1 - hosts[i, j+1] ) ) + w1[j] * hosts[i, j+1] # calculate individual fitnesses
            end
        end

        for i = 1:L-1   # calculate parasite fitnesses
            wp0[i] = 1 - s*host1[i]
            wp1[i] = 1 - s*host0[i]
        end

        # host reproduction and selection
        wMax = maximum(w)
        for i = 1:N
            w[i] = w[i] / wMax  # calculate relative fitnesses
        end
        
        for i = 1:N     # select random host haplotypes to go to next generation
            randomNum = rand(1:N)
            if(rand() <= w[randomNum])
                hostsnew[i, :] = hosts[randomNum, :]
            end
        end
        global hosts = hostsnew
        # mutation
        for i = 1:N
            for j = 1:L # loop through loci in all individuals
                if(rand() <= u) # if random num is less than mutation rate
                    hosts[i, j] = abs(hosts[i, j] - 1)  # flip the allele of the locus
                end
            end
        end
        # recombination
        for i = 1:Int(N/2)
            if(hosts[2*i, 1] == 1 || hosts[2*i - 1, 1] == 1) # loop through pairs of hosts, checking for recombination allele
                for j = 1:L-1
                    if(rand() <= r) # 0.5 chance at each locus to swap between pairs
                        global temp = hosts[2*i, :]
                        hosts[2*i, j+1:L] = hosts[2*i - 1, j+1:L]
                        hosts[2*i - 1, j+1:L] = temp[j+1:L]
                    end
                end
            end
        end
        
        # parasites
        for i = 1:L-1   # calculate parasite fitness
            wp0[i] = 1 - s*host1[i+1]
            wp1[i] = 1 - s*host0[i+1]
        end
        for pgen = 1:npgens # guess who finally figured out dot operater thing
            wp_bar .= wp0 .* para0 .+ wp1 .* para1  # mean parasite fitness
            para0 .= para0 .* wp0 ./ wp_bar         # freq. of parasite wild type alleles
            para1 .= 1 .- para0                     # freq. of parasite mutant alleles

            para0 .= para0 .* (1 - u) .+ para1 .* u # wild type allele freq. after mutation
            para1 .= 1 .- para0                     # mutant type allele freq. after mutation
        end

        # fixation
        for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
            tempSum = 0
            for j = 1:N
                tempSum += hosts[j, i]
            end
            host1[i] = tempSum/N
            host0[i] = 1 - host1[i]
        end
        hostMutantAlleleFreq[rep, gen, :] .= host1  # store mutant allele freq.
        frecomb[rep] = host1[1]     # store ending freq. of recombination allele
        
        if(!fix && host1[1] > 0.99)
            global fix = true
            global nfix += 1
            gfix[rep] = gen
            break
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
        genrecomb[gen] += host1[1]/nreps         # add this generation's recomb allele freq. to genrecomb
        genrecombarr[gen, rep] = host1[1]        # add this generation's recomb allele freq. to genrecombarr
    end
end

genrecombstd = Array{Float64}(undef, ngens)
for i = 1:ngens
    genrecombstd[i] = std(genrecombarr[i, :])
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
@show genrecombarr
@show genrecombstd
@show maximum(genrecombarr)

using Plots
testPlot = plot(1:ngens, genrecomb, yerror = genrecombstd)
display(testPlot)
readline()