#=
This program is for fitting the dynode curve to extract the fluorescence time
Every file has the following structure: [time], [voltage], [empty]
=#

## Import libraries
begin
    using LsqFit
    using CSV, DataFrames
    using Unitful, UnitfulRecipes
    using StatsPlots
    using LaTeXStrings
    using Measurements
    using LinearAlgebra
    using NaNStatistics # MOVING MEAN
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
    df[!, :Height] = abs.(df[!, :Height] * u"V") .|> u"mV"
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


## Getting some errors for the data
# plot data with errors
# plot(df[!, :Time], df[!, :Height] .± df[!, :σ], label="Data with errors")
# plot!(df[!, :Time], movmean(df[!, :Height], length(df[!, :Height]) / 50), label="Moving average")


## LsqFit     Leo 1994  7.2 -> pg 158
strt = 370;
t = df[!, :Time][strt:end] .- df[!, :Time][strt];
h = df[!, :Height][strt:end];
w = df[!, :σ][strt:end].^(-2);
N(t, p) = @. p[1]/p[2] * exp(-t/p[2]);*
p0 = [120.0, 230.0];

# fit = curve_fit(N, ustrip.(t), ustrip.(h), ustrip.(w), p0); ## with errors
fit = curve_fit(N, ustrip.(t), ustrip.(h), p0); ## without errors
params = fit.param;
errrors = stderror(fit)
fitted_function(t) = N(ustrip.(t), params);


## Plot with before after comparison
save_loc = "plots\\NaI(Ti)_ResponseTime_Fluoresence_Decay.pdf"

p_ez = begin
    plot(ustrip.(t .+ df[!, :Time][strt]), h, color=:blue, label=:none);
    plot!(df[!, :Time][begin:strt], df[!, :Height][begin:strt], color=:blue, label="Collected data"); # plot data

    annotate!(275, 20, text(L"N = \frac{N_0}{\tau_d} \exp\left(\frac{-t}{\tau_d}\right)", :red, 10))
    plot!(t .+ df[!, :Time][strt], N(ustrip.(t), ustrip.(params)),
            color=:red, label= "Fit to data, decay time: $((params[2] ± errrors[2])*u"ns")",
            # xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7
    ) # add fit

    # With error:
    # every = 20
    # scatter!(df[!, :Time][1:every:end], (df[!, :Height] .± df[!, :σ])[1:every:end],
    #         label="Selected errors", markersize=3, color=:magenta, dpi=500
    # )
    # savefig(save_loc)
end

## Advanced fit -> Leo 1994  7.2 -> pg 158
N_7_3(t, p) = @. p[1] * exp(-t/p[2]) + p[3] * exp(-t/p[4])

p0 = [500.0, 230.0, 500.0, 230.0];
lb = zeros(length(p0));

fit2 = curve_fit(N_7_3, ustrip.(t), ustrip.(h), ustrip(w), p0, lower=lb);
adv_params = fit2.param
adv_errrors = stderror(fit2) # LAPACK exception

fitted_function(t) = N_7_3(ustrip.(t), params)

p_adv = begin
    plot(t .+ df[!, :Time][strt], h, color=:blue, label=:none);
    plot!(df[!, :Time][begin:strt], df[!, :Height][begin:strt], color=:blue, label="Collected data"); # plot data

    annotate!(375, 15, text(L"N = A \cdot \exp\left(\frac{-t}{\tau_f}\right) \! + \! B \cdot \exp\left(\frac{-t}{\tau_s}\right)", :red, 10))
    plot!(t .+ df[!, :Time][strt], N_7_3(ustrip.(t), ustrip.(adv_params)),
            color=:red, label= "Fit to data, decay time: $((adv_params[2] ± 8.9)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7,
            axisfontsize=4
    ) # add fit
    # scatter!(df[!, :Time][1:every:end], (df[!, :Height] .± df[!, :σ])[1:every:end], label="Selected errors", markersize=3)
    # savefig(save_loc)
end

## plot both parts
#
# l = @layout [a{0.5h};
#              b{0.5h}]

# plot(p_ez, p_adv, layout=l) # dpi=900)

##
pmtt = df[!, :Time];
pmth = df[!, :Height];
pmte = df[!, :σ];
pmtw = (1 ./ pmte).^2;
N_pmt(t, p) = @. p[1] * exp(-t/p[2]) + p[3]*exp(-t/p[4]) - p[5] * exp(-t/p[6])


p_pmt = begin
    p0pmt = [200.0, 230.0, 80.0, 5.0, 50, 5.0];
    # lbpmt = zeros(length(p0pmt))
    lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2];

    # rise_start = 328   ## This is the standard one

    rise_start = 366 # no error -> 252
    # rise_start = 366 # w error

    # fit_pmt = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), ustrip.(pmtw[rise_start:end]), p0pmt, lower=lbpmt)
    fit_pmt = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), p0pmt, lower=lbpmt)
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
            legendfontsize=7,
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

    scatter!([pmtt[rise_start]], [pmth[rise_start]], color=:darkred, label="Plot from which i start the red fit")
    # annotate!(750, 75, text(L"""A e^{-\frac{t}{\tau_s}} + B e^{-\frac{t}{\tau_f}} - C e^{-\frac{t}{\tau}}""", :red, 10));
    # annotate!(750, 60, text(latstring[1], 9) );
    # annotate!(750, 50, text(latstring[2], 9) );
    # annotate!(750, 40, text(latstring[3], 9) );
    annotate!(325, 40, text(L"""A e^{-\frac{t}{\tau_s}} + B e^{-\frac{t}{\tau_f}} - C e^{-\frac{t}{\tau}}""", :red, 10));
    annotate!(325, 23, text(latstring[1], 9) );
    annotate!(325, 13, text(latstring[2], 9) );
    annotate!(325, 3, text(latstring[3], 9) );
end

begin
    l = @layout [a{0.35h};
                b{0.65h}]
    plot(p_ez, p_pmt, layout=l, dpi=1000)
end
# savefig("plots\\3.1Fluorescence_7_3&8.13.pdf")



using Bootstrap

function bootstrap_func(d1, d2, i)
    tmp = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), ustrip.(pmtw[rise_start:end]), p0pmt, lower=lbpmt)
    return tmp.param[i]
end


rise_start = 340
tb = ustrip.(pmtt[rise_start:end])
th = ustrip.(pmth[rise_start:end])
bsdata = (tb, th)
n_boot = 1000

# bs1 = bootstrap(bootstrap_func, bsdata, BalancedSampling(n_boot))

##
# function jacobian_model(x,p)
#         J = Array{Float64}(undef, length(x), length(p))
#         @. J[:,1] = exp(-x*p[2])
#         @. @views J[:,2] = x * p[1] / (p[2]^2)*J[:,1]
#
#         @. J[:,3] = exp(-x*p[4])
#         @. @views J[:,4] = x * p[3] / (p[4]^2)*J[:,3]
#
#         @. J[:,5] = exp(-x*p[6])
#         @. @views J[:,6] = x * p[5] / (p[5]^2) * J[:,5]
#
#         J
# end
# # model_inplace(F, x, p) = (@. F = N_pmt(x, p));
#
# fit = curve_fit(N_pmt, jacobian_model,
#                 ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]),
#                 p0pmt;
#                 lower=lbpmt
# )
#
# for i in 1:1:450
#     testfit = fit = curve_fit(N_pmt, jacobian_model,
#                     ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]),
#                     p0pmt;
#                     lower=lbpmt
#     )
#     try
#         printstyled(i , ": "; bold=true, color=:cyan)
#         test_covar = stderror(testfit)
#         printstyled(string(test_covar); bold=true, color=:red)
#         print("\n")
#     catch e
#         println(e)
#     end
# end
