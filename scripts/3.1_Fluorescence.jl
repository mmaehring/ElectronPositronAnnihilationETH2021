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

## Loading data
begin
    path = "data\\3.1\\Dynode3.1\\800V\\800V_ALL0000_CurveData.txt"
    df = CSV.read(path, DataFrame; select=[1,2]);
    # unproc_t = df[!, :Time]
    # unproc_v = df[!, :Height]
    df[!, :Time] = df[!, :Time] * 1e9; # nanoseconds
    df[!, :Time] = df[!, :Time] .- minimum(df[!, :Time]);
    df[!, :Time] = df[!, :Time] * u"ns"
    df[!, :Height] = abs.(df[!, :Height] * u"V") .|> u"mV" # * -1
    # df[!, :Height] = -(df[!, :Height] * u"V") .|> u"mV" # * -1
    df[!, :σ] .= abs.(df[!, :Height] .- 0.25*movmean(df[!, :Height], length(df[!, :Height]) / 50)) / 10;
    # df[!, :σ] .= 0.05 * df[!, :Height] .+ 0.0001u"mV"
end;



# Initial Plot and moving average to extract average of error
begin
    plot(df[!, :Time], df[!, :Height], size=(900, 600), dpi=900)
    i = 50
    smoothed_data = movmean(df[!, :Height], length(df[!, :Height]) / i)
    plot!(df[!,:Time], smoothed_data, label="Smoothing with window of $(i)")
end

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


## Getting some errors for the data
# plot data with errors
# plot(df[!, :Time], df[!, :Height] .± df[!, :σ], label="Data with errors")
# plot!(df[!, :Time], movmean(df[!, :Height], length(df[!, :Height]) / 50), label="Moving average")


## LsqFit     Leo 1994  7.2 -> pg 158
strt = 370;
t = df[!, :Time][strt:end] .- df[!, :Time][strt];
h = df[!, :Height][strt:end];
w = df[!, :σ][strt:end].^(-2);

N(t, p) = @. p[1]/p[2] * exp(-t/p[2]);
N_shift(t, p) = @. p[1]/p[2] * exp(-(t-p[3])/p[2]);

p0 = [120.0, 230.0];

# fit = curve_fit(N, ustrip.(t), ustrip.(h), ustrip.(w), p0); ## with errors
fit = curve_fit(N, ustrip.(t), ustrip.(h), p0) ## without errors
params = fit.param;
errrors = stderror(fit);
fitted_function(t) = N(ustrip.(t), params);

## Plot with before after comparison
save_loc = "plots\\NaI(Ti)_ResponseTime_Fluoresence_Decay.pdf"

p_ez = begin
    plot(ustrip.(t .+ df[!, :Time][strt]), h, color=:blue, label=:none);
    plot!(df[!, :Time][begin:strt], df[!, :Height][begin:strt], color=:blue, label="Collected data"); # plot data

    annotate!(275, 20, text(L"N = \frac{N_0}{\tau_d} \exp\left(\frac{-t}{\tau_d}\right)", :red, 10))
    plot!(t .+ df[!, :Time][strt], N(ustrip.(t), ustrip.(params)),
            color=:red, label= "Fit to data, decay time: $((params[2] ± errrors[2])*u"ns")",
            ylabel="Height (mV)", margin=3mm,
            legendfontsize=7
    )
    # With error:
    # every = 20
    # scatter!(df[!, :Time][1:every:end], (df[!, :Height] .± df[!, :σ])[1:every:end],
    #         label="Selected errors", markersize=3, color=:magenta, dpi=500
    # )
    # savefig(save_loc)
end

## PMT inclusion
begin
    pmtt = df[!, :Time];
    pmth = df[!, :Height];
    pmte = df[!, :σ];
    pmtw = (1 ./ pmte).^2;
end
N_pmt(t, p) = @. p[1] * exp(-t/p[2]) + p[3]*exp(-t/p[4]) - p[5] * exp(-t/p[6])
# N_pmt(t, p) = @. p[1] * exp((t - p[7])/p[2]) + p[3]*exp((t-p[7])/p[4]) - p[5] * exp((t-p[7])/p[6]) # shift the time with a bias

p_pmt = begin
    p0pmt = [200.0, 230.0, 80.0, 5.0, 50, 5.0];
    # p0pmt = [200.0, 230.0, 80.0, 5.0, 50, 5.0, 104.0]; # w shift
    # lbpmt = zeros(length(p0pmt))

    lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0.0];
    # lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2];
    # ubpmt = [Inf, 500.0, Inf, 500.0, Inf, 500.0, 200]; # for the bias

    rise_start = 328   ## This is the standard one

    # rise_start = 366 # no error -> 252
    # rise_start = 366 # w error
    # rise_start = 194 # w shift

    # fit_pmt = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), ustrip.(pmtw[rise_start:end]), p0pmt, lower=lbpmt)
    fit_pmt = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), p0pmt, lower=lbpmt)#, upper=ubpmt)
    pmt_params = fit_pmt.param
    # pmt_error = stderror(fit_pmt)
    pmt_errors = :none

    # if pmt_error != :none
    #     @show pmt_error
    # end
    # @show round.(pmt_params; digits=2)

    pmt_fit = plot(pmtt, pmth, color=:blue, label="Collected data");
    plot!(pmtt[1:end], N_pmt(ustrip.(pmtt[1:end]), Measurements.value.(ustrip.(pmt_params))),
            color=:red, label= "Fit to data",#, decay time: $(round((pmt_params[1]), digits=2)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7, xguidefontsize=9, yguidefontsize=9,
            axisfontsize=3, dpi=800, ylims=(-1,125)
    ) #
    latstring = let p = round.(pmt_params; digits=2)
        A = p[1]
        τs = p[2]
        B = p[3]
        τf = p[4]
        C = p[5]
        τ = p[6]
        (L"""\quad A = %$(A), \ \tau_s = %$(τs)""", L"""B = %$(B), \ \ \  \tau_f = %$(τf)""", L"""C=%$(C), \ \ \tau=%$(τ)""")
    end

    # scatter!([pmtt[rise_start]], [pmth[rise_start]], color=:darkred, label="Plot from which i start the red fit")
    # annotate!(750, 75, text(L"""A e^{-\frac{t}{\tau_s}} + B e^{-\frac{t}{\tau_f}} - C e^{-\frac{t}{\tau}}""", :red, 10));
    # annotate!(750, 60, text(latstring[1], 9) );
    # annotate!(750, 50, text(latstring[2], 9) );
    # annotate!(750, 40, text(latstring[3], 9) );
    annotate!(325, 40, text(L"""A e^{-\frac{t}{\tau_s}} + B e^{-\frac{t}{\tau_f}} - C e^{-\frac{t}{\tau}}""", :red, 10));
    annotate!(325, 23, text(latstring[1], 9) );
    annotate!(325, 13, text(latstring[2], 9) );
    annotate!(325, 3, text(latstring[3], 9) );
    # annotate!(750, 80, "$(round(pmt_params[end], digits=2))" );
    # savefig("plots\\")
end

begin
    l = @layout [a{0.35h};
                b{0.65h}]
    plot(p_ez, p_pmt, layout=l, dpi=1000)
    savefig("plots\\3.1Fluorescence_7_1&7_3&8.13.pdf")
end

## Bootstrapping
const perform_bootstrap = false
if perform_bootstrap
    nboot= 1
    N_pmt(t, p) = @. p[1] * exp(-t/p[2]) + p[3]*exp(-t/p[4]) - p[5] * exp(-t/p[6])
    xdata = ustrip.(pmtt[rise_start:end])
    ydata = ustrip.(pmth[rise_start:end])

    bs_params = zeros((nboot, 6))
    lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2];
    begin
        for i in ProgressBar(1:nboot)
            p0pmt = [200.0 + 5*randn(), 130.0+5*rand(), 80.0 + 0.5*randn(), 5.0 + 0.1*randn(), 50 + 0.5*randn(), 5.0 + 0.1*randn()];
            # p0pmt = [200.0, 230.0, 80.0, 5.0, 50, 5.0];

            sample = StatsBase.sample(1:length(xdata), length(xdata); replace=true)
            bsdata = (xdata[sample], ydata[sample])

            fit_bs = curve_fit(N_pmt, bsdata[1], bsdata[2], p0pmt, lower=lbpmt)
            bs_param = fit_bs.param
            bs_params[i,1:end] = bs_param
            @info(bs_param)
        end
        writedlm("processed_data\\Bootstrap_$(nboot)_values_PMT_Fit_Starting_$(rise_start).txt", bs_params)
    end
end

# File to get bootstrap
path_366_bs = "processed_data\\Bootstrap_100_000_values_PMT_Fit_NotStarting0_366.txt"
bootstrapped_366_df = CSV.read(path_366_bs, DataFrame)
bs_error366 = std(bootstrapped_366_df |> Array, dims=1);
pmt_params_366 = [400.64540996899206, 252.99482289817277, 79.99934922173671,  4.757826783723315, 819.8722666253608,  84.71862037780365];

path_328_bs = "processed_data\\Bootstrap_100_000_values_PMT_Fit_Starting_328.txt"
bootstrapped_328_df = CSV.read(path_328_bs, DataFrame)
bs_error328 = std(bootstrapped_328_df |> Array, dims=1);
pmt_params_328 = [274.3454791596396, 298.3392776658812, 79.99967595242731, 4.110641790077652, 7074.019077968472, 36.65530625785385];
# @info 100*diag(errors_of_the_mean_params ./ pmt_params_328)


function evaluate_fit_and_bootstrap(bs_errors, xdata, ydata, σdata, params, fct)
    printstyled("Starting evaluation \n", bold=true, color=:red)

    @info "Fitted params:"
    @info params
    @info bs_errors

    @info "Relative errors of the fitted parameters:"
    @info 100 * diag(bs_errors ./ params)


    fctN(x) = fct(x, params)
    @info "χ² characteristics in the following order: "
    @info "dof, χ², χ²dof"
    @info χ²_calc(fctN, ustrip.(xdata), ustrip.(ydata), ustrip.(σdata), length(params))


    # @info
    printstyled("End of information", bold=true, color=:red)
end

evaluate_fit_and_bootstrap(bs_error328, pmtt[328:end], pmth[328:end], pmte[328:end], pmt_params_328, N_pmt)
# evaluate_fit_and_bootstrap(bs_error366, pmtt[366:end], pmth[366:end], pmte[366:end], pmt_params_366, N_pmt)


## Fit to smoothed data -> didnt work
# begin
#     path = "data\\3.1\\Dynode3.1\\800V\\800V_ALL0000_CurveData.txt"
#     df = CSV.read(path, DataFrame; select=[1,2]);
#     # unproc_t = df[!, :Time]
#     # unproc_v = df[!, :Height]
#
#     df[!, :Time] = df[!, :Time] * 1e9; # nanoseconds
#     df[!, :Time] = df[!, :Time] .- minimum(df[!, :Time]);
#     df[!, :Time] = df[!, :Time] * u"ns"
#     df[!, :Height] = abs.(df[!, :Height] * u"V") .|> u"mV"
#     df[!, :σ] .= abs.(df[!, :Height] .- 0.25*movmean(df[!, :Height], length(df[!, :Height]) / 50)) / 10;
#     # df[!, :σ] .= 0.05 * df[!, :Height] .+ 0.0001u"mV"
# end;
#
#
# # Initial Plot and moving average to extract average of error
# begin
#     plot(df[!, :Time], df[!, :Height], size=(900, 600), dpi=900)
#     i = 20
#     smoothed_data = movmean(df[!, :Height], length(df[!, :Height]) / i)
#     plot!(df[!,:Time], smoothed_data, label="Smoothing with window of $(i)")
# end
#
# N_pmt_smod(t, p) = @. p[1] * exp((t-p[7])/p[2]) + p[3]*exp((t-p[7])/p[4]) - p[5] * exp((t-p[7])/p[6])
# p0smod = [200.0, 230.0, 80.0, 5.0, 50, 5.0, 110.0];
# # p0pmt = [200.0, 230.0, 80.0, 5.0, 50, 5.0, 104.0]; # w shift
# # lbpmt = zeros(length(p0pmt))
#
# lbsmod = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0.0, 1e-2];
# ubsmod = [5000, 500.0, Inf, 5000.0, Inf, 5000.0, 200]; # for the bias
#
# fit_smod = curve_fit(N_pmt_smod, ustrip.(df[!, :Time]), ustrip.(smoothed_data), p0smod, lower=lbsmod, upper=ubsmod)
# smod_params = fit_smod.param
#
# begin
#     plot(df[!, :Time], df[!, :Height], size=(900, 600), dpi=900, margin=3mm)
#     plot!(df[!,:Time], smoothed_data, label="Smoothing with window of $(i)")
#     plot!(df[!; :Time], )
# end
