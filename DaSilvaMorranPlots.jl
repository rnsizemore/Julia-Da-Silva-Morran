using CSV
using DataFrames
using Plots
# call syntax
# julia .\DaSilvaMorranPlots.jl 2025 01 26 0.1 1 2 1.0 1000 30 3 2 1
# julia .\DaSilvaMorranPlots.jl year month day pnd recomb treatment h nreps ngens nscens fitscen
year = ARGS[1]
month = ARGS[2]
day = ARGS[3]
pnd = parse(Float64,ARGS[4])
recomb = parse(Int,ARGS[5])
treatment = parse(Int,ARGS[6])
h = parse(Float64,ARGS[7])
nreps = parse(Int,ARGS[8])
ngens = parse(Int,ARGS[9])
nscens = parse(Int,ARGS[10])
num = parse(Int, ARGS[11])
#pnd, recomb, treatment, h, nreps, ngens = 0.1, 1, 2, 1.0, 100, 30

# check if file name will have count
if num != 1
    recombcsvname = "$(year)-$(month)-$(day)_genrecomb_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    recombcsvpath = joinpath(@__DIR__,  "$(recombcsvname).csv")
    outcrosscsvname = "$(year)-$(month)-$(day)_outcross_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    outcrosscsvpath = joinpath(@__DIR__,  "$(outcrosscsvname).csv")
    malefitcsvname = "$(year)-$(month)-$(day)_malefit_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    malefitcsvpath = joinpath(@__DIR__,  "$(malefitcsvname).csv")
    hermfitcsvname = "$(year)-$(month)-$(day)_hermfit_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    hermfitcsvpath = joinpath(@__DIR__,  "$(hermfitcsvname).csv")
else
    recombcsvname = "$(year)-$(month)-$(day)_$(num)_genrecomb_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    recombcsvpath = joinpath(@__DIR__,  "$(recombcsvname).csv")
    outcrosscsvname = "$(year)-$(month)-$(day)_$(num)_outcross_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    outcrosscsvpath = joinpath(@__DIR__,  "$(outcrosscsvname).csv")
    malefitcsvname = "$(year)-$(month)-$(day)_$(num)_malefit_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    malefitcsvpath = joinpath(@__DIR__,  "$(malefitcsvname).csv")
    hermfitcsvname = "$(year)-$(month)-$(day)_$(num)_hermfit_pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)_ngens$(ngens)_nreps$(nreps)";
    hermfitcsvpath = joinpath(@__DIR__,  "$(hermfitcsvname).csv")
end

# bring in data as a matrix
recomb_data = CSV.File(recombcsvpath) |> Tables.matrix
outcross_data = CSV.File(outcrosscsvpath) |> Tables.matrix
malefit_data = CSV.File(malefitcsvpath) |> Tables.matrix
hermfit_data = CSV.File(hermfitcsvpath) |> Tables.matrix

# initalize averages array (3 rows for the 3 scenarios)
recomb_averages = zeros(Float64, ngens+1, 3)
outcross_averages = zeros(Float64, ngens+1, 3)
malefit_averages = zeros(Float64, ngens, 3)
hermfit_averages = zeros(Float64, ngens, 3)
for i = 1:nscens
    # save every third row of data to get only scen i
    recomb_scen_data = recomb_data[i:nscens:(nreps*nscens), 3:(ngens+3)]
    outcross_scen_data = outcross_data[i:nscens:(nreps*nscens), 3:(ngens+3)]
    malefit_scen_data = malefit_data[i:nscens:(nreps*nscens), 3:(ngens+2)]
    hermfit_scen_data = hermfit_data[i:nscens:(nreps*nscens), 3:(ngens+2)]
    for j = 1:nreps
        recomb_averages[:, i] .= recomb_averages[:, i] .+ recomb_scen_data[j, :]
        outcross_averages[:, i] .= outcross_averages[:, i] .+ outcross_scen_data[j, :]
        malefit_averages[:, i] .= malefit_averages[:, i] .+ malefit_scen_data[j, :]
        hermfit_averages[:, i] .= hermfit_averages[:, i] .+ hermfit_scen_data[j, :]
    end
end
# divide values of averages array by nreps to calculate average
recomb_averages .= recomb_averages ./ nreps
outcross_averages .= outcross_averages ./ nreps
malefit_averages .= malefit_averages ./ nreps
hermfit_averages .= hermfit_averages ./ nreps

dir_name = "pnd$(pnd)_recomb$(recomb)_treatment$(treatment)_h$(h)"
if !isdir(dir_name)
    mkpath(dir_name)
end

gr()
labels = reshape(map((x,y)->"N=$x, r=$y",[1500, 750, 750],[0.5, 0.5, 0.1]),1,nscens)

# outcrossing plot
outcrossPlot = plot(0:ngens, outcross_averages,label=labels,title="Outcrossing Rate, $(nreps) Replicates",
                titlefont = font(14), xlabel="Generation",ylabel="Outcrossing Frequency",ylims=(0,1));
#display(outcrossPlot)
outcross_plot_file = "outcross_ngens$(ngens)_nreps$(nreps)"
savefig(joinpath(dir_name, "$(outcross_plot_file).pdf"))
# recombination allele plot
recombPlot = plot(0:ngens, recomb_averages,label=labels, title="Recombination Allele Frequency, $(nreps) Replicates",
                titlefont = font(14), xlabel="Generation",ylabel="Recombination Allele Frequency",ylims=(0,1));
recomb_plot_file = "genrecomb_ngens$(ngens)_nreps$(nreps)"
savefig(joinpath(dir_name, "$(recomb_plot_file).pdf"))
#display(recombPlot)
# male and hermaphrodite fitness plot
# fitplot = nothing # clear existing plot
fitscen = parse(Int, ARGS[12])
fitplot = plot(1:ngens, malefit_averages[:, fitscen],label="Males",title="Male and Hermaphrodite Fitness, $(nreps) Replicates",
                titlefont = font(14), xlabel="Generation",ylabel="Average Fitness",ylims=(0,1));
plot!(hermfit_averages[:, fitscen], label = "Hermaphrodites")
fit_plot_file = "fitness_ngens$(ngens)_nreps$(nreps)"
savefig(joinpath(dir_name, "$(fit_plot_file).pdf"))
#display(fitplot)