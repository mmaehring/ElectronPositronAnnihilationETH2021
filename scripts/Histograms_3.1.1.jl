using StatsBase
using CSV
using DataFrames
using Measurements
using StatsPlots
using StatsBase
using Unitful
using UnitfulRecipes
using Measures

path ="/home/marcus/Desktop/PraktikumIII/e+e-A/3.3.1_TMP/"
cd(path)
searchdir(path,key) = filter(x->contains(x,key), readdir(path))
file_names = searchdir(path, ".txt")

begin
        mv300 = file_names[1]
        mv200I = file_names[end-3]
        mv200E = file_names[end-2]

        # mv300 = CSV.read(mv300, DataFrame)
        mv200I = CSV.read(mv200I, DataFrame)
        mv200E = CSV.read(mv200E, DataFrame)

        # mv300[!, :σ] = sqrt.(mv300[!,:Count])
        mv200I[!, :σ] = sqrt.(mv200I[!,:Count])
        mv200E[!, :σ] = sqrt.(mv200E[!,:Count])
end

begin
        p1 = plot(mv200I[!, "LLD"], mv200I[!, "Count"], label=:none, xlims=(-0.075,2.6), xlabel="LLD Voltage (V)", ylabel="Counts");
        scatter!(mv200I[!, "LLD"], mv200I[!, "Count"] .± mv200I[!, "σ"],
                label="200mV, 850V, I", color="red", markersize=2.5);

        p2 = plot(mv200E[!, "LLD"], mv200E[!, "Count"], label=:none, xlims=(-0.075,2.6), xlabel="LLD Voltage (V)");
        scatter!(mv200E[!, "LLD"], mv200E[!, "Count"] .± mv200E[!, "σ"],
                        label="50mV, 850V, E", color="red", markersize=2.5);

        l = @layout [a{0.5w} b{0.5w}]
        plot(p1,p2, layout=l, dpi=800, size = (1200, 500), legendfontsize=10,
                margin=5mm
        )
        # savefig("Scattered_Points_3.3.1.png")
end
## "Histograms"

begin
        h1 = plot(mv200I[!, "LLD"], mv200I[!, "Count"], seriestype=:steppost, xlabel="LLD Voltage (V)", ylabel="Counts",
                  xlims=(-0.05,2.6), label="Counts per interval, in [0.01, 0.5]V",
                  xticks=0:0.2:2.6,
        );
        scatter!(mv200I[!, "LLD"].+0.025, mv200I[!, "Count"].± mv200I[!, "σ"],
                markersize=1.25, color="blue", label=:none,
                title="850V bias, 200mV window, 10 s acquisition time"
        );

        # for the histogram we want the point to be in the center of the step
        # mv200E[!, "LLD"] = mv200E[!, "LLD"] .+ 0.05/2
        h2 = plot(mv200E[!, "LLD"], mv200E[!, "Count"], seriestype=:steppost, xlabel="LLD Voltage (V)",
                  xlims=(-0.05,2.55), label="Counts per interval: (50±3)mV",
                  xticks=0:0.25:2.50
        );
        scatter!(mv200E[!, "LLD"].+0.025, mv200E[!, "Count"].± mv200E[!, "σ"],
                markersize=1.25, color="blue", label=:none,
                title="850V bias, 50mV window, 20 s acquisition time"
        );
        l = @layout [a{0.5w} b{0.5w}];
        plot(h1, h2, layout=l, dpi=800, size = (1400, 500),
             legend=:topright, legendfontsize=10,
             margin = 5mm,
             titlefontsize=12
        )
        # savefig("Histogram_3.3.1.png")
end
