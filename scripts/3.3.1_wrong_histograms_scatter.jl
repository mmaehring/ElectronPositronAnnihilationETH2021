using DataFrames
using Plots
using CSV
using Statistics
using Unitful
using StatsPlots
using StatsBase

## Imprecise ζ
ζvoltages_LLD_1 = [0  , 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.45,
                0.5, 0.6, 0.7, 0.8,  1, 1.5, 2, 2.5, 3.0, 3.5, 4.0]
ζcounts_1 = [795, 741, 651, 403, 333, 253, 177, 108,  209, 423,  561,
           504, 178,  56,  44, 51,  19, 2,   1,   0,   0,   0]

ζvoltages_LLD_2 = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,
                   0.31, 0.32, 0.35, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.5, 0.51,
                   0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.62, 0.65, 0.7, 0.8, 0.9, 1]

ζcounts_2 = [112, 114, 143, 108,  14,  16, 15,  13, 14,  10, 9, 3, 3, 3, 3, 2, 2,
            3, 0, 2, 11, 13, 10, 19, 20, 24, 41, 35, 21, 16, 20, 15, 2, 1, 1, 0, 0]

plot(ζvoltages_LLD_1, ζcounts_1, linestyle=:dash);
xticks_vals = collect.([0:0.25:2, 2:0.5:4]);
xticks_vals = vcat(xticks_vals[1], xticks_vals[2]);
scatter!(ζvoltages_LLD_1, ζcounts_1, color="red",
            pointstyle=".", title="Inprecise Window: 200 mV, 10s measuring - 1",
            xlabel="LLD Threshhold (V)",
            ylabel="Detected counts",
            pointsize=0.25,
            xticks=xticks_vals,
            yticks=0:100:850,
            titlefont=12
)
# savefig("200mV10s1.png")

plot(ζvoltages_LLD_2, ζcounts_2, linestyle=:dash);
xticks_vals = collect.([0:0.25:2, 2:0.5:4]);
xticks_vals = vcat(xticks_vals[1], xticks_vals[2]);
scatter!(ζvoltages_LLD_2, ζcounts_2, color="red",
            pointstyle=".", title="Inprecise: Window=200 mV, 10s measuring - 2",
            xlabel="LLD Threshhold (V)",
            ylabel="Detected counts",
            pointsize=0.25,
            xticks=xticks_vals,
            yticks=0:25:175,
            titlefont=12
)
# savefig("200mV10s2.png")

## Precise ξ
cd("C:/Users/marcu/OneDrive/Desktop/PraktikumIII/e+e-_Annihilation/data/3.3/3.3.1/")
searchdir(path,key) = filter(x->contains(x,key), readdir(path));
runs = searchdir("PreciseMeasurements", ".txt");

ξ300mV = CSV.read("PreciseMeasurements"*"/"*runs[1], DataFrame)
ξ600mV = CSV.read("PreciseMeasurements"*"/"*runs[2], DataFrame)

plot(ξ300mV[!, :LLD ], ξ300mV[!, :Count], linestyle=:dash);
xticks_vals = collect.([0:0.25:2, 2:5]);
xticks_vals = vcat(xticks_vals[1], xticks_vals[2]);
scatter!(ξ300mV[!, :LLD ], ξ300mV[!, :Count], color="red",
        pointstyle=".", title="Precise: Window=300 mV, 20s measuring",
        xlabel="LLD Threshhold (V)",
        ylabel="Detected counts",
        pointsize=0.25,
        xticks=xticks_vals,
        yticks=0:250:2250,
        titlefont=12
)
# savefig("300mV20s.png")

plot(ξ600mV[!, :LLD ], ξ600mV[!, :Count], linestyle=:dash);
xticks_vals = collect.([0:0.5:2, 2:5]);
xticks_vals = vcat(xticks_vals[1], xticks_vals[2]);
scatter!(ξ600mV[!, :LLD ], ξ600mV[!, :Count], color="red",
        pointstyle=".", title="Precise: Window=600 mV, 20s measuring",
        xlabel="LLD Threshhold (V)",
        ylabel="Detected counts",
        pointsize=0.5,
        xticks=xticks_vals,
        yticks=0:250:3000,
        titlefont=12
)
# savefig("600mV20s.png")

## With histograms:



# 200mV10 - 2
histogram(ζvoltages_LLD_2, bins=0:0.2:1.0, weights=ζcounts_2,
        title="Inprecise: Window=200 mV, 10s measuring - 2",
        xlabel="LLD Threshhold (V)",
        ylabel="Detected counts",
        pointsize=0.25,
        titlefont=12,
        xticks=0:0.2:1.0
)
# savefig("200mV10s_hist.png")


## 300mV20
# h = fit(Histogram, ξ300mV[!, :LLD], 0:0.3:1.5)

histogram(ξ300mV[!, :LLD], bins=0:0.3:1.5, weights=ξ300mV[!, :Count],
        title="Precise: Window=300 mV, 20s measuring",
        xlabel="LLD Threshhold (V)",
        ylabel="Detected counts",
        titlefont=12,
        xticks=0:0.3:1.5
)
# savefig("300mV20s_hist.png")
a = [ξ300mV[!, :Count][findall(x->i<x<=i+0.3, ξ300mV[!, :LLD])]  for i in [0, 0.3, 0.6, 0.9, 1.2]]
weight = mean.(a) # better counts?
