#=
This program is for plotting the different anode curves from part 3.1. to compare
for saturation.
Data names given in form AnodeData[voltage]V_[persist time]s.txt (comma separated)
Every file has the following structure: [time], [voltage], [empty]
The empty was generated from the oscilloscope and we don't need to remove it so here it stays :)
=#

## Loading packages
begin
    using CSV, DataFrames
    using Unitful, UnitfulRecipes
    using StatsPlots
    using LaTeXStrings
end


## Loading data
begin
    ## Loading 800 V 5 s
    path800V = "data\\3.1\\Anode3.1\\800V\\ALL0003\\AnodeData800V_INFs.txt"
    D800V_df = CSV.read(path800V, DataFrame, select=[1,2]);
    D800V_df[!, :time] = (-D800V_df[!, :time][1] .+ D800V_df[!, :time]) * 1e9; # nanoseconds

    ## Loading 900 V 5 s
    path900V = "data\\3.1\\Anode3.1\\AnodeSaturation_V1\\900V_5s\\AnodeData900V_5s.txt"
    D900V_df = CSV.read(path900V, DataFrame, select=[1,2]);
    D900V_df[!, :time] = (-D900V_df[!, :time][1] .+ D900V_df[!, :time]) * 1e9; # nanoseconds

    ## Loading 910 V 5 s
    path910V = "data\\3.1\\Anode3.1\\AnodeSaturation_V1\\910V_5s\\AnodeData910V_5s.txt"
    D910V_df = CSV.read(path910V, DataFrame, select=[1,2]);
    D910V_df[!, :time] = (-D910V_df[!, :time][1] .+ D910V_df[!, :time]) * 1e9; # nanoseconds

    ## Loading 1000 V
    path1000V = "data\\3.1\\Anode3.1\\AnodeSaturation_V1\\1000V_5s\\AnodeData1000V_5s.txt";
    D1000V_df = CSV.read(path1000V, DataFrame, select=[1,2]);
    D1000V_df[!, :time] = (-D1000V_df[!, :time][1] .+ D1000V_df[!, :time]) * 1e9; # nanoseconds
end;


## Trying to make a big mess plotting everything at the same time

begin
    plot(D800V_df[!, :time], D800V_df[!, :voltage], label="800V")
    plot!(D900V_df[!, :time], D900V_df[!, :voltage], label="900V")
    plot!(D910V_df[!, :time], D910V_df[!, :voltage], label="910V")
    plot!(D1000V_df[!, :time], D1000V_df[!, :voltage], label="1000V")
end

#=
Tragically, this did not work. Even though the oscilloscope images the data captured
the last 5 seconds (or since the start) it only saves the last measurement, meaning the
saturation doesn't occur in the data (unless you are ridiculously lucky) :/
=#
