include("DaSilvaMorranDistributed.jl")

# order of parameters
# julia .\DaSilvaMorranJob.jl workers nreps treatment s pnd ngens pnd h h_recomb recomb_scen

# base parameters:
# julia .\DaSilvaMorranJob.jl 5 100 2 0.8 0.0001 30 0.1 1 1 1
#                      workers^       s^   u^          h^   ^recomb_scen
#                          nreps^           ngens^        ^h_recomb
#                          treatment^             pnd^        
@everywhere using Plots
@everywhere using Dates
@everywhere using Statistics
@everywhere using CSV
@everywhere using DataFrames
@everywhere using Random
@everywhere using Distributions
@everywhere using StatsBase

# define nreps from user input
nreps = parse(Int,ARGS[2])

# parameters to test
treatment = parse(Int,ARGS[3])       # 0 - control 1 - evolution 2 - coevolution
s = parse(Float64,ARGS[4])             # max selection coeff.
u = parse(Float64,ARGS[5])           # mutation rate per locus
ngens = parse(Int,ARGS[6])          # num. of generations
pnd = parse(Float64,ARGS[7])           # prob. of x-chromosome nondisjunction
h = parse(Float64,ARGS[8])               # dominance of alleles
h_recomb = parse(Float64,ARGS[9])      # dominance of recombination allele
recomb_scen = parse(Int,ARGS[10])     # 0 - recomb. allele has no effect 1 - recomb. allele has effect
nscens = 3

@everywhere function work_function(treatment, nrep, s, u, ngens, pnd, h, h_recomb, recomb_scen, nscens)
    IL = 30             # num. of host interaction loci
    L = IL + 1          # num. of host loci
    initfr = 0.2        # inital freq. of males
    standing = true     # create standing genetic variation
    NArray = [1500, 750, 750]            # host pop.
    rArray = [0.5, 0.5, 0.1]             # recombination rate
    npgens = 1          # parasite gen. per host gen.

    # variables
    hosts = zeros(Int, NArray[1], L)                            # holds value for each locus in each genome
    hostsnew = zeros(Float64, NArray[1], L)                         # holds hosts value for next generation
    phenarr0 = zeros(Float64, NArray[1], L)                         # holds value for phenotype of each locus in each individual, wild type 0
    phenarr1 = zeros(Float64, NArray[1], L)                         # holds value for phenotype of each locus in each individual, wild type 1
    temp = Array{Int}(undef, L)                                 # temp holder for an individuals genes
    host0 = zeros(L)                                            # freq. of mutant allele in population at each locus
    host1 = zeros(L)                                            # freq. of wild type allele in population at each locus
    phen0 = zeros(L)
    phen1 = zeros(L)
    para0 = zeros(L-1)                                          # freq. of parasite mutant allele
    para1 = zeros(L-1)                                          # freq. of parasite wild type allele
    w0 = Array{Float64}(undef, L)                               # host allele fitness for mutant alleles
    w1 = Array{Float64}(undef, L)                               # host allele fitness for wild type alleles
    wp0 = Array{Float64}(undef, L-1)                            # parasite mutant allele fitnesses
    wp1 = Array{Float64}(undef, L-1)                            # parasite wild type allele fitnesses
    wp_bar = Array{Float64}(undef, L-1)                         # parasite mean fitness
    w_male = zeros(Float64, ngens, nscens)
    w_herm = zeros(Float64, ngens, nscens)
    fix = false                                                 # flag for fixation of recombination allele
    nfix = 0                                                    # number of simulations where recombination allele fixates
    gfix = zeros(Int)                                    # stores the generation where recombination allele comes to fixation
    frecomb = Array{Float64}(undef)                      # freq. of recombination allele at end of simulation run
    indiv_vector = zeros(Int, 750)                              # vector of individual sexes
    genrecomb = zeros(Float64, ngens+1, nscens)                 # average freq. of recombination allele at each generation over the replicates
    genrecombarr = zeros(Float64, ngens+1, nscens)       # stores the array of values of recomb. freq. at each generation
    genrecombstd = zeros(Float64, ngens+1, nscens)                    # stores standard deviation of replicates recomb. freq.
    malefreq = zeros(Float64, ngens+1, nscens)                  # average male freq. at each generation over the replicates
    malefreqarr = zeros(Float64, ngens+1, nscens)       # stores the array of values of male. freq. at each generation
    malefreqstd = zeros(Float64, ngens+1, nscens)                    # stores standard deviation of replicates male freq.
    outcross = zeros(Float64, ngens+1, nscens)                  # average outcrossing rate at each generation over the replicates
    outcrossarr = zeros(Float64, ngens+1, nscens)       # stores the array of values of outcrossing rate at each generation
    outcrossstd = zeros(Float64, ngens+1, nscens)                    # stores standard deviation of replicates outcrossing rate
    hostMutantAlleleFreq = Array{Float64}(undef, ngens+1, L)   # host mutant allele freq. for each generation of a simulation run
    meanf = Array{Float64}(undef, ngens+1, L)                         # mean host mutant allele freq. over all simulation runs
    iFixCount = zeros(Int, L-1)                          # number of simulation runs where each locus has fixation
    iPrevFix = Array{Int}(undef, L)                             # the last fixed allele for each host interaction locus
    gen = 1             # current generation count 
    scen = 1            # current scenario
    pgen = 1            # current parasite generation count
    #random_genomes = rand(Float64, Int(N[1]/2))    # random number for inital freq of mutant allele

    for scen in 1:3
        # local initialization steps for each scenario
        N = NArray[scen]
        r = rArray[scen]
        print("Scenario $scen\n")
        print("N=$N, r=$r\n")
        nfix = 0
        gfix = zeros(Int)
        hosts = zeros(Int, N, L)
        hostsnew = zeros(Int, N, L)
        phenarr0 = zeros(Float64, N, L)
        phenarr1 = zeros(Float64, N, L)

        n_indivs = round(Int,N/2)
        w = Array{Float64}(undef, n_indivs)
        n_males = round(Int, ceil(n_indivs * initfr))
        indiv_vector = zeros(Int, n_indivs)
        indiv_vector[1:n_males] .= 1 # set inital fraction of males to the same frequency of recomb allele 
        shuffle!(indiv_vector)
        freqm = sum(indiv_vector)/n_indivs
        hosts[:, 1] =  rand(Bernoulli(initfr),N)
        if standing
            hosts[:, 2:L] = rand(0:1, N, IL)
        else
            hosts[:, 2:L] .= 0 # wildtype allele state is 0
        end
        for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
            host1[i] = sum(hosts[:, i])/N
            host0[i] = 1.0 - host1[i]
        end
        # initialization of data collection structures
        gen = 1
        genrecomb[gen, scen] += host1[1]          # add generation 0's recomb allele freq. to genrecomb
        malefreq[gen, scen] += freqm # add generation 0's male freq. to malefreq
        outcross[gen, scen] += 2*(freqm-pnd) # add generation 0 outcrossing rate
        #if genrecomb[gen,scen]>1 || genrecomb[gen,scen]<0
        #    error("Out of range value in genrecomb[$gen,$scen]")
        #end
        genrecombarr[gen, scen] = host1[1]        # add generation 0's recomb allele freq. to genrecombarr
        hostMutantAlleleFreq[gen, :] = host1  # store mutant allele freq. at generation 0
        # intitalize wild type allele at fixation in parasite population
        para0 .= 1
        para1 .= 0

        fix = false                 # reset fixation of recomb allele
        iFixCount = zeros(Int, L)
        iPrevFix = zeros(Int, L)  # intitalize last fixed allele as 0

        for gen = 1:ngens # loop through generations
            #print("\n Gen #", gen, "\n")
            #print("Start of gen: ", sum(hosts[:, 1])/N, "\n")
            if treatment != 0
                @. w0[2:L] = 1 - s*para0  # wild type allele fitness
                @. w1[2:L] = 1 - s*para1  # mutant allele fitness
                # generate phenotype arr
                for i = 1:n_indivs
                    ind_i = i * 2 - 1
                    for j = 2:L
                        phenarr1[i, j] = 2*h*middle(hosts[ind_i, j], hosts[ind_i+1, j]) + (1-2*h)*hosts[ind_i, j]*hosts[ind_i+1, j]
                    end
                    for j = 2:L
                        phenarr0[i, j] = 2*h*middle((1-hosts[ind_i, j]), (1-hosts[ind_i+1, j])) + (1-2*h)*(1-hosts[ind_i, j])*(1-hosts[ind_i+1, j])
                    end
                    # generate phenotype for recomb allele
                    phenarr1[i, 1] = 2*h_recomb*middle(hosts[ind_i, 1], hosts[ind_i+1, 1]) + (1-2*h_recomb)*hosts[ind_i, 1]*hosts[ind_i+1, 1]
                    phenarr0[i, 1] = 2*h_recomb*middle((1-hosts[ind_i, 1]), (1-hosts[ind_i+1, 1])) + (1-2*h_recomb)*(1-hosts[ind_i, 1])*(1-hosts[ind_i+1, 1])
                end
                for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
                    phen1[i] = sum(phenarr1[:, i])/n_indivs
                    phen0[i] = 1.0 - phen1[i]
                end
                #phen0 .= sum(phenarr0[:, 2:L])/n_indivs
                #phen1 .= sum(phenarr1[:, 2:L])/n_indivs
                #print("start fitness calc")
                w = ones(n_indivs)          # intitalize individual fitnesses
                for i = 1:n_indivs
                    #@show phenarr1
                    #@show phenarr2
                    for j = 2:L
                        w[i] = w[i] * ( w0[j] * phenarr0[i, j] + w1[j] * phenarr1[i, j] ) # calculate individual fitnesses
                    end
                end
                #print("finish fitness calc")
            else
                w = ones(N)
            end
            #@show maximum(w)
            # calculate parasite fitnesses
            @. wp0 = 1.0 - s*phen1[2:L]
            @. wp1 = 1.0 - s*phen0[2:L]

            # host reproduction and selection
            wMax = maximum(w)
            @. w = w / wMax   # calculate relative fitnesses
            #print("started loop\n")
            #@show indiv_vector
            temp_vector = zeros(Int, n_indivs)
            hermaphrodite_ind = findall(x -> x == 0, indiv_vector)
            male_ind = findall(x -> x == 1, indiv_vector)
            w_male[gen, scen] += mean(w[male_ind])
            w_herm[gen, scen] += mean(w[hermaphrodite_ind])
            for i = 1:2:N
                ind_i = round(Int,ceil(i / 2))
                # index of sampled hermaphrodite
                #print("start h sample")
                h_sample_ind = sample(hermaphrodite_ind, Weights(w[hermaphrodite_ind]))
                #print("finish h sample")
                # index of sampled hermaphrodite genome
                h_sample = 2 * (h_sample_ind - 1) + rand(1:2)
                #print("finish h sample\n")
                if rand() < (freqm)       # freqm chance for outcrossing, otherwise selfing
                    #print("start outcrossing\n")
                    # index of sampled male
                    m_sample_ind = sample(male_ind, Weights(w[male_ind]))
                    # index of sampled male genome
                    m_sample = 2 * (m_sample_ind - 1) + rand(1:2) 
                    #print("finish m sample\n")
                    # outcrossing between sampled genomes
                    if(recomb_scen == 1 && (hosts[h_sample, 1] == 1 || hosts[m_sample, 1] == 1) || recomb_scen == 0)
                        hostsnew[i, 1] = hosts[h_sample, 1]
                        hostsnew[i+1, 1] = hosts[m_sample, 1]
                        for j = 1:L-1
                            if(rand() <= r) # 0.5 chance at each locus to swap between pairs
                                hostsnew[i, j+1:L] = hosts[m_sample, j+1:L]
                                hostsnew[i + 1, j+1:L] = hosts[h_sample, j+1:L]
                            end
                        end
                    else
                        hostsnew[i, :] = hosts[h_sample, :]
                        hostsnew[i+1, :] = hosts[m_sample, :]
                    end
                    if rand() <= (0.5 + pnd)
                        temp_vector[ind_i] = 1
                    else
                        temp_vector[ind_i] = 0
                    end
                    #print("finish outcrossing\n")
                else
                    #selfing
                    # recombination between 2 genomes of individual
                    #print("start selfing\n")
                    h_sample_2 = 2 * (h_sample_ind - 1) + rand(1:2) # sampling of 2nd hermaphrodite genome
                    if h_sample == h_sample_2
                        hostsnew[i, :] = hosts[h_sample, :]
                        hostsnew[i+1, :] = hosts[h_sample, :]
                    else 
                        if(recomb_scen == 1 && (hosts[h_sample, 1] == 1 || hosts[h_sample_2, 1] == 1) || recomb_scen == 0)
                            hostsnew[i, 1] = hosts[h_sample, 1] # taking values of recomb. allele from parents
                            hostsnew[i+1, 1] = hosts[h_sample_2, 1]
                            for j = 1:L-1
                                if(rand() <= r) # 0.5 chance at each locus to swap between pairs
                                    hostsnew[i, j+1:L] = hosts[h_sample_2, j+1:L]
                                    hostsnew[i + 1, j+1:L] = hosts[h_sample, j+1:L]
                                end
                            end
                        else
                            hostsnew[i, :] = hosts[h_sample, :]
                            hostsnew[i+1, :] = hosts[h_sample_2, :]
                        end
                    end
                    if(rand() <= pnd)   # pnd chance for male
                        temp_vector[ind_i] = 1
                    else
                        temp_vector[ind_i] = 0
                    end
                    #print("finish selfing\n")
                end
                #print("completed ", (i+1)/2, " loops\n")
            end
            indiv_vector = temp_vector
            hosts = hostsnew
            freqm = sum(indiv_vector)/n_indivs
            #print("finished loop\n")
            #print("\n Wild type fitness: ", w0[1], "\n")
            #print(" Wild type parasite fitness: ", wp0[1], "\n \n")
            #for i = 1:N
            #    while(true)
            #        j = rand(1:N)
            #        if(rand() <= w[j])
            #            hostsnew[i, :] = hosts[j, :]
            #            break
            #        end
            #    end
            #end
            
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
            #for i = 1:2:N
            #    if(hosts[i, 1] == 1 || hosts[i + 1, 1] == 1) # loop through pairs of hosts, checking for recombination allele
            #        for j = 1:L-1
            #            if(rand() <= r) # 0.5 chance at each locus to swap between pairs
            #                temp = hosts[i, :]
            #                hosts[i, j+1:L] = hosts[i+1, j+1:L]
            #                hosts[i + 1, j+1:L] = temp[j+1:L]
            #            end
            #        end
            #    end
            #end
            #print("After recombination: ", sum(hosts[:, 1])/N, "\n")
            # parasites
            if treatment == 2
                for pgen = 1:npgens
                    #for i = 1:L-1
                        @. wp_bar = wp0 * para0 + wp1 * para1  # mean parasite fitness
                        @. para0 = para0 * wp0 / wp_bar         # freq. of parasite wild type alleles
                        @. para1 = 1.0 - para0                     # freq. of parasite mutant alleles

                        @. para0 = para0 * (1.0 - u) + para1 * u # wild type allele freq. after mutation
                        @. para1 = 1.0 - para0                     # mutant type allele freq. after mutation
                    #end
                end
            end

            # fixation
            for i = 1:L     # each item in host1 is freq. of mutant allele at each locus, each item host0 is wild type allele freq. at each locus
                host1[i] = sum(hosts[:, i])/N
                host0[i] = 1.0 - host1[i]
            end
            hostMutantAlleleFreq[gen+1, :] = host1  # store mutant allele freq.
            frecomb = host1[1]     # store ending freq. of recombination allele
            
            if(!fix && host1[1] > 0.99)
                fix = true
                nfix += 1
                gfix = gen
                #break
            end
            
            for i = 2:L
                if(host1[i] > 0.99 && iPrevFix[i] == 0)
                    iFixCount[:] = iFixCount[:] .+ 1
                    iPrevFix[i] = 1
                end
                if(host0[i] > 0.99 && iPrevFix[i] == 1)
                    iFixCount[:] = iFixCount[:] .+ 1
                    iPrevFix[i] = 0
                end
            end
            #print("before: ", genrecomb[gen+1, scen], "\n")

            outcross[gen+1, scen] += 2*(freqm-pnd)
            malefreq[gen+1, scen] += freqm
            genrecomb[gen+1, scen] += host1[1]         # add this generation's recomb allele phenotype freq. to genrecomb
            #if genrecomb[gen+1,scen]>1 || genrecomb[gen+1,scen]<-1
            #    error("Out of range value in genrecomb[$(gen+1),$scen]")
            #end
            outcrossarr[gen+1, scen] = 2*(freqm-pnd)
            genrecombarr[gen+1, scen] = host1[1]        # add this generation's recomb allele phenotype freq. to genrecombarr
            malefreqarr[gen+1, scen] = freqm
            #print("End of gen: ", sum(hosts[:, 1])/N, "\n")
        end

        for i = 0:ngens
            outcrossstd[i+1] = std(outcrossarr[i+1, scen])
            genrecombstd[i+1] = std(genrecombarr[i+1, scen])
            malefreqstd[i+1] = std(malefreqarr[i+1, scen])
        end
        println("\nFixations at Host Interaction Loci")
        println("No. of fixations/locus/gen: ",  sum(iFixCount) / (IL) / (ngens+1))
        println("Apparent no. of fixations/locus/gen: ", count(hostMutantAlleleFreq[ngens+1, 2:L] .> 0.99) / (IL) / (ngens+1), "\n")

        println("Fixation of Recombination Allele")
        println("Prop. of replicates with fixation: ", nfix)
        if(nfix > 0)
            println("Mean no. of generations to fixation: ", sum(gfix) / nfix)
        end
        println("Mean freq of recomb allele: ", sum(frecomb))
        println("\n")
    end
    # return variables we want to save
    return (nrep=nrep, genrecomb1=genrecomb[:, 1], genrecomb2=genrecomb[:, 2], genrecomb3=genrecomb[:, 3], outcross1=outcross[:, 1], outcross2=outcross[:, 2], outcross3=outcross[:, 3], male1=w_male[:, 1], male2=w_male[:, 2], male3=w_male[:, 3], herm1=w_herm[:, 1], herm2=w_herm[:, 2], herm3=w_herm[:, 3])
end

# start running the event loop on workers
start_workers()

# array of job calls
jobs = []
for i = 1:nreps
    push!(jobs, (:work_function, treatment, i, s, u, ngens, pnd, h, h_recomb, recomb_scen, nscens))
end

submit_jobs(jobs)

n = length(jobs)

date = Dates.format(now(), "YYYY-mm-dd-HH-MM")
plotname = "pnd$(pnd)_recomb$(recomb_scen)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)"

recombcsvfile = "$(date)_genrecomb_$(plotname).csv"
outcrosscsvfile = "$(date)_outcross_$(plotname).csv"
malefitcsvfile = "$(date)_malefit_$(plotname).csv"
hermfitcsvfile = "$(date)_hermfit_$(plotname).csv"

# generate column labels
colnames = Any["scen", "nrep"]
for gen = 0:ngens
    push!(colnames, "gen$(gen)")
end
# fitness data has no generation 0, so remove that column name
fitcolnames = deleteat!(copy(colnames), 3)

for i in 1:n
    results = take!(resultsqueue)
    @info "Obtained result ($i/$n)" pairs(results)
    for scen in 1:nscens
        CSV.write(recombcsvfile, DataFrame(permutedims(vcat([scen,results.nrep],results[1+scen])), colnames); append=isfile(recombcsvfile))
        CSV.write(outcrosscsvfile, DataFrame(permutedims(vcat([scen,results.nrep],results[4+scen])), colnames); append=isfile(outcrosscsvfile))
        CSV.write(malefitcsvfile, DataFrame(permutedims(vcat([scen,results.nrep],results[7+scen])), fitcolnames); append=isfile(malefitcsvfile))
        CSV.write(hermfitcsvfile, DataFrame(permutedims(vcat([scen,results.nrep],results[10+scen])), fitcolnames); append=isfile(hermfitcsvfile))
    end
end
stop_workers()
#=
# reading csv
# genrecomb = Matrix(CSV.read("nreps100_ngens30_control.CSV", DataFrame))
dir_name = "pnd$(pnd)_recomb$(recomb_scen)_treatment$(treatment)_h$(h)"
if !isdir(dir_name)
    mkpath(dir_name)
end

labels = reshape(map((x,y)->"N=$x, r=$y",NArray,rArray),1,nscens)
gr()
# outcrossing plot
outcrossPlot = plot(0:ngens, outcross,label=labels,plot_title="$nreps reps per scenario, pnd = $pnd",
                xlabel="Generation",ylabel="Outcrossing Frequency",ylims=(0,1));
#display(outcrossPlot)
outcross_plot_file = "outcross_ngens$(ngens)_nreps$(nreps)"
savefig(joinpath(dir_name, "$(outcross_plot_file).pdf"))
# recombination allele plot
recombPlot = plot(0:ngens, genrecomb,label=labels,plot_title="$nreps reps per scenario, pnd = $pnd",
                xlabel="Generation",ylabel="Average Fitness",ylims=(0,1));
recomb_plot_file = "genrecomb_ngens$(ngens)_nreps$(nreps)"
savefig(joinpath(dir_name, "$(recomb_plot_file).pdf"))
#display(recombPlot)
# male and hermaphrodite fitness plot
testscen = 1    # scenario to plot
fitplot = plot(1:ngens, w_male[:, testscen],label=labels,plot_title="$nreps reps per scenario, pnd = $pnd",
                xlabel="Generation",ylabel="Average Fitness",ylims=(0,1));
plot!(w_herm[:, testscen])
fit_plot_file = "fitness_ngens$(ngens)_nreps$(nreps)_scen$(testscen)"
savefig(joinpath(dir_name, "$(fit_plot_file).pdf"))
#display(fitplot)


#csv_file = "nreps$(nreps)_ngens$(ngens)_$(treatment)_indiv_test"
#CSV.write("$(csv_file).CSV", Tables.table(genrecomb))
=#