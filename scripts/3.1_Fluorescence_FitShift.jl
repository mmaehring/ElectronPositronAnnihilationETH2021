#=
This program is for fitting the dynode curve to extract the fluorescence time
Every file has the following structure: [time], [voltage], [empty]
=#

## Import libraries
begin
    using Logging
    using LsqFit
    using CSV, DataFrames
    using Unitful, UnitfulRecipes, StatsPlots, LaTeXStrings, Measures
    using Measurements
    using LinearAlgebra
    using Statistics
    import StatsBase
    using ProgressBars, DelimitedFiles
end

cd(raw"C:\Users\marcu\OneDrive\Desktop\PraktikumIII\e+e-_Annihilation")

## Helper function

function χ²_calc(fit_function, x, y, σ, dependent_vars = 2)
    if σ == :none
        χ² = sum((fit_function.(x) .- y).^2)
        dof = length(x) - dependent_vars - 1
        χ²dof = χ²/dof
        return dof, χ², χ²dof
    end
    χ² = sum((fit_function.(x) .- y).^2 ./ σ.^2)
    dof = length(x) - dependent_vars - 1
    χ²dof = χ²/dof
    return dof, χ², χ²dof
end

## Loading data
begin
    ### Raw data
    path = "data\\3.1\\Dynode3.1\\800V\\800V_ALL0000_CurveData.txt"
    df = CSV.read(path, DataFrame; select=[1, 2]);
    df[!, :Time] = df[!, :Time] * 1e9; # nanoseconds
    df[!, :Time] = df[!, :Time] * u"ns"
    # df[!, :Time] = df[!, :Time] .+ 370u"ns"
    df[!, :Height] = (df[!, :Height] * u"V") .|> u"mV"
    df[!, :Height] = df[!, :Height] * -1
    df[!, :σ] .= 0.05*1u"mV"
end;



p0pmt = [274.0, 298.0, 80.0, 4.11, 7074, 36.66];
lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0.0];
NPMT(t, p) = @. p[1] * exp((t-400)/p[2]) + p[3]*exp((t-400)/p[4]) - p[5] * exp((t-400)/p[6])


begin
    plot(df[!, :Time], df[!, :Height], legend=:bottomright, label="Data")
    vline!([0], label="Start of signal")
    fitted_data = NPMT(df[!, :Time].|>ustrip, p0pmt)
    plot(df[!, :Time], fitted_data)
end


# NPMT(t, p) = @. p[1] * exp((t-p[7])/p[2]) + p[3]*exp((t-p[7])/p[4]) - p[5] * exp((t-p[7])/p[6])
fit_pmt = curve_fit(NPMT, ustrip.(pmtt), ustrip.(pmth), p0pmt, lower=lbpmt)



begin
    pmt_fit = plot(pmtt, pmth, color=:blue, label="Collected data");
    plot!(pmtt[1:end], N_pmt(ustrip.(pmtt[1:end]), Measurements.value.(ustrip.(pmt_params))),
    color=:red, label= "Fit to data",#, decay time: $(round((pmt_params[1]), digits=2)*u"ns")",
    xlabel="Time series",
    ylabel="Height (mV)",
    legendfontsize=7, xguidefontsize=9, yguidefontsize=9,
    axisfontsize=3, dpi=800, ylims=(-1,125)
    ) #
    # latstring = let p = round.(pmt_params; digits=2)
    #     A = p[1]
    #     τs = p[2]
    #     B = p[3]
    #     τf = p[4]
    #     C = p[5]
    #     τ = p[6]
    #     (L"""\quad A = %$(A), \ \tau_s = %$(τs)""", L"""B = %$(B), \ \ \  \tau_f = %$(τf)""", L"""C=%$(C), \ \ \tau=%$(τ)""")
    # end
    #
    # annotate!(325, 40, text(L"""A e^{-\frac{t}{\tau_s}} + B e^{-\frac{t}{\tau_f}} - C e^{-\frac{t}{\tau}}""", :red, 10));
    # annotate!(325, 23, text(latstring[1], 9) );
    # annotate!(325, 13, text(latstring[2], 9) );
    # annotate!(325, 3, text(latstring[3], 9) );
    # savefig("plots\\")
end
