using CSV
using DataFrames
using Plots
using Statistics
# call syntax
# julia .\DaSilvaMorranPlots.jl 2025 02 11 21 42 0.1 1 2 1.0 1000 30 3 2 1 2
# julia .\DaSilvaMorranPlots.jl year month day hour minute pnd recomb treatment h nreps ngens nscens fitscen model
year = ARGS[1]
month = ARGS[2]
day = ARGS[3]
hour = ARGS[4]
minute = ARGS[5]
pnd = parse(Float64,ARGS[6])
recomb = parse(Int,ARGS[7])
treatment = parse(Int,ARGS[8])
h = parse(Float64,ARGS[9])
nreps = parse(Int,ARGS[10])
ngens = parse(Int,ARGS[11])
nscens = parse(Int,ARGS[12])
model = parse(Int, ARGS[14]) # 1 - new model 2 - old model
#pnd, recomb, treatment, h, nreps, ngens = 0.1, 1, 2, 1.0, 100, 30

date = "$(year)-$(month)-$(day)-$(hour)-$(minute)"

if model == 1
    plotname = "pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)"

    recombcsvname = "$(date)_genrecomb_$(plotname)";
    recombcsvpath = joinpath(@__DIR__,  "$(recombcsvname).csv")
    outcrosscsvname = "$(date)_outcross_$(plotname)";
    outcrosscsvpath = joinpath(@__DIR__,  "$(outcrosscsvname).csv")
    malefitcsvname = "$(date)_malefit_$(plotname)";
    malefitcsvpath = joinpath(@__DIR__,  "$(malefitcsvname).csv")
    hermfitcsvname = "$(date)_hermfit_$(plotname)";
    hermfitcsvpath = joinpath(@__DIR__,  "$(hermfitcsvname).csv")

    # bring in data as a DataFrame
    recomb_data = CSV.read(recombcsvpath, DataFrame)
    outcross_data = CSV.read(outcrosscsvpath, DataFrame)
    malefit_data = CSV.read(malefitcsvpath, DataFrame)
    hermfit_data = CSV.read(hermfitcsvpath, DataFrame)
    recomb_averages = DataFrame()
    outcross_averages = DataFrame()
    malefit_averages = DataFrame()
    hermfit_averages = DataFrame()
else
    plotname = "treatment$(treatment)_ngens$(ngens)_nreps$(nreps)"

    recombcsvname = "$(date)_old_genrecomb_$(plotname)";
    recombcsvpath = joinpath(@__DIR__,  "$(recombcsvname).csv")
    
    recomb_data = CSV.read(recombcsvpath, DataFrame)
    recomb_averages = DataFrame()
end

col_names = collect(Symbol("gen$t") for t in 0:ngens)
recomb_averages = combine(groupby(recomb_data,:scen), col_names .=> mean .=> col_names)
if model == 1
    outcross_averages[:, col] = combine(groupby(outcross_data,:scen), col_names .=> mean .=> col_names)
    # exclude gen 0 for fitness
    fit_col_names = collect(Symbol("gen$t") for t in 1:ngens)
    malefit_averages[:, col] = combine(groupby(malefit_data,:scen), fit_col_names .=> mean .=> fit_col_names)
    hermfit_averages[:, col] = combine(groupby(hermfit_data,:scen), fit_col_names .=> mean .=> fit_col_names)
end

dir_name = "pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)"
if !isdir(dir_name)
    mkpath(dir_name)
end

gr()
labels = reshape(map((x,y)->"N=$x, r=$y",[1500, 750, 750],[0.5, 0.5, 0.1]),1,nscens)

# plots does not support dataframes by default, so its converted to a Matrix
# potential todo: use StatsPlots package, which does support dataframes
# also potential todo: change how averages dataframe is made to avoid having to use permutedims()

# recombination allele plot
recombPlot_data = Matrix(permutedims(recomb_averages[:, 2:32]))
recombPlot = plot(0:ngens,recombPlot_data,label=labels, title="Recombination Allele Frequency, $(nreps) Replicates",
                titlefont = font(14), xlabel="Generation",ylabel="Recombination Allele Frequency",ylims=(0,1));
recomb_plot_file = "genrecomb_ngens$(ngens)_nreps$(nreps)"
savefig(joinpath(dir_name, "$(recomb_plot_file).pdf"))

if model == 1
    # outcrossing plot
    outcrossPlot_data = Matrix(permutedims(outcross_averages[:, 2:32]))
    outcrossPlot = plot(0:ngens,outcrossPlot_data,label=labels,title="Outcrossing Rate, $(nreps) Replicates",
                    titlefont = font(14), xlabel="Generation",ylabel="Outcrossing Frequency",ylims=(0,1));
    #display(outcrossPlot)
    outcross_plot_file = "outcross_ngens$(ngens)_nreps$(nreps)"
    savefig(joinpath(dir_name, "$(outcross_plot_file).pdf"))

    # male and hermaphrodite fitness plot
    # fitscen is which scenario to plot
    fitscen = parse(Int, ARGS[13])
    malefitplot_data = Matrix(permutedims(malefit_averages[:, 2:31]))
    hermfitplot_data = Matrix(permutedims(hermfit_averages[:, 2:31]))
    fitplot = plot(1:ngens,malefitplot_data[:, fitscen],label="Males",title="Male and Hermaphrodite Fitness, $(nreps) Replicates",
                    titlefont = font(14), xlabel="Generation",ylabel="Average Fitness",ylims=(0,1));
    plot!(hermfitplot_data[:, fitscen], label = "Hermaphrodites")
    fit_plot_file = "fitness_ngens$(ngens)_nreps$(nreps)"
    savefig(joinpath(dir_name, "$(fit_plot_file).pdf"))
    #display(fitplot)
end