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
    using Statistics, NaNStatistics # MOVING MEAN
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
    df[!, :Time] = df[!, :Time] .+ 370u"ns"
    df[!, :Height] = (df[!, :Height] * u"V") .|> u"mV"
    df[!, :Height] = df[!, :Height] * -1
    df[!, :σ] .= 0.05*1u"mV"
end;

begin
    plot(df[!, :Time], df[!, :Height], legend=:bottomright, label="Data")
    vline!([0], label="Start of signal")
end

## LsqFit     Leo 1994  7.2 -> pg 158

complete_fit = begin
    p0pmt = [200.0, 230.0, 80.0, 5.0, 50, 5.0];
    lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0];


    rise_start = 151
    pmtt = df[!, :Time][rise_start:end]
    pmth = df[!, :Height][rise_start:end]

    tₛₕ = +9
    N_pmt(t, p) = @. p[1] * exp(-(t+tₛₕ)/p[2]) + p[3]*exp(-(t+tₛₕ)/p[4]) - p[5] * exp(-(t+tₛₕ)/p[6])
    fit_pmt = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), p0pmt, lower=lbpmt)
    pmt_params = fit_pmt.param
    try
        pmt_error = stderror(fit_pmt)
        println(pmt_error)
    catch
        println(":(")
    end

    pmt_fit = plot(df[!,:Time], df[!, :Height], color=:blue, label="Collected data");
    # scatter!((df[!,:Time][rise_start].|>ustrip, 0))
    # scatter!(-15,0)
    plot!(pmtt[rise_start+4:end], N_pmt(ustrip.(pmtt[rise_start+4:end]), Measurements.value.(ustrip.(pmt_params))),
            color=:red, label= "Fit",#, decay time: $(round((pmt_params[1]), digits=2)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7, xguidefontsize=9, yguidefontsize=9,
            axisfontsize=3.5, dpi=800, ylims=(-1,125)
    ) #
    latstring = let p = round.(pmt_params; digits=2)
        A = p[1]
        τs = p[2]
        B = p[3]
        τf = p[4]
        C = p[5]
        τ = p[6]
        (L"""\quad A = %$(A), \ \tau_s = %$(τs)""", L"""B = %$(B), \ \ \  \tau_f = %$(τf)""", L"""C=%$(C), \ \ \tau=%$(τ) """, L"""t_s=%$(tₛₕ) """)
    end

    annotate!(630, 48+20, text(L"""A e^{-\frac{t-t_s}{\tau_s}} + B e^{-\frac{t-t_s}{\tau_f}} - C e^{-\frac{t-t_s}{\tau}}""", :red, 12));
    annotate!(180, 43, text(latstring[1], 9) );
    annotate!(180, 33, text(latstring[2], 9) );
    annotate!(180, 23, text(latstring[3], 9) );
    annotate!(180, 13, text(latstring[4], 9) );
    # savefig("plots\\3.1PMT_ShiftedFit.pdf")
end

simple_fit = begin
    strt = 425;
    t = df[!, :Time][strt:end];
    h = df[!, :Height][strt:end];
    w = df[!, :σ][strt:end].^(-2);

    p0 = [120.0, 230.0];
    t_s = 65

    N(t, p) = @. p[1]/p[2] * exp(-(t+t_s)/p[2]);


    fit = curve_fit(N, ustrip.(t), ustrip.(h), p0) ## without errors
    params = fit.param;
    errrors = stderror(fit);
    fitted_function(t) = N(ustrip.(t), params);

    # Plotting
    plot(df[!, :Time], df[!, :Height], color=:blue, label="Collected data"); # plot data
    annotate!(230, 20, text(L"N = \frac{N_0}{\tau_d} \exp\left(\frac{-(t-t_s)}{\tau_d}\right)", :red, 10))
    plot!(t .+ df[!, :Time][strt], N(ustrip.(t), ustrip.(params)),
            color=:red, label= "Fit, decay time: $((params[2] ± errrors[2])*u"ns")",
            ylabel="Height (mV)", margin=3mm, xguidefontsize=9, yguidefontsize=9,
            legendfontsize=7
    )
end

begin
    l = @layout [a{0.35h};
                b{0.65h}]
    plot(simple_fit, complete_fit, layout=l, dpi=1000)
    # savefig("plots\\3.1Fluorescence_7_1&7_3&8.13_SHIFTED.pdf")
end
