## Import libraries
begin
    using LsqFit
    using CSV, DataFrames
    using Unitful, UnitfulRecipes
    using StatsPlots
    using LaTeXStrings
    using Measurements
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
end;

# plot(unproc_t*1e9, unproc_v*1e3)
# vline!([-360])

## Initial Plot
strt = 370;
t = df[!, :Time][strt:end] .- df[!, :Time][strt];
h = df[!, :Height][strt:end];
plot(t, h)

# strt = 573;
# t = df[!, :Time][strt:end] .- df[!, :Time][strt];
# h = df[!, :Height][strt:end];
# plot(t, h)

## LsqFit     Leo 1994  7.2 -> pg 158
N(t, p) = @. p[1]/p[2] * exp(-t/p[2]);
p0 = [120.0, 230.0];

fit = curve_fit(N, ustrip.(t), ustrip.(h), p0);
params = fit.param;
errrors = stderror(fit);
fitted_function(t) = N(ustrip.(t), params);


## Plot with before after comparison
save_loc = "plots\\NaI(Ti)_ResponseTime_Fluoresence_Decay.pdf"

p_ez = begin
    plot(ustrip.(t .+ df[!, :Time][strt]), h, color=:blue, label=:none);
    plot!(df[!, :Time][begin:strt], df[!, :Height][begin:strt], color=:blue, label="Collected data"); # plot data

    annotate!(275, 15, text(L"N = \frac{N_0}{\tau_d} \exp\left(\frac{-t}{\tau_d}\right)", :red, 10))
    plot!(t .+ df[!, :Time][strt], N(ustrip.(t), ustrip.(params)),
            color=:red, label= "Fit to data, decay time: $((params[2] ± 6.5)*u"ns")",
            # xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7
    ) # add fit
    # savefig(save_loc)
end

## Advanced fit -> Leo 1994  7.2 -> pg 158
N_adv(t, p) = @. p[1] * exp(-t/p[2]) + p[3] * exp(-t/p[4])

p0 = [500.0, 230.0, 500.0, 230.0];
lb = zeros(length(p0));

fit2 = curve_fit(N_adv, ustrip.(t), ustrip.(h), p0, lower=lb);
adv_params = fit2.param
adv_errrors = confidence_interval(fit2) # LAPACK exception

fitted_function(t) = N_adv(ustrip.(t), params)

p_adv = begin
    plot(t .+ df[!, :Time][strt], h, color=:blue, label=:none);
    plot!(df[!, :Time][begin:strt], df[!, :Height][begin:strt], color=:blue, label="Collected data"); # plot data

    annotate!(375, 15, text(L"N = A \cdot \exp\left(\frac{-t}{\tau_f}\right) \! + \! B \cdot \exp\left(\frac{-t}{\tau_s}\right)", :red, 10))
    plot!(t .+ df[!, :Time][strt], N_adv(ustrip.(t), ustrip.(adv_params)),
            color=:red, label= "Fit to data, decay time: $((adv_params[2] ± 8.9)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7,
            axisfontsize=4
    ) # add fit
    # savefig(save_loc)
end

## plot both parts

l = @layout [a{0.5h};
             b{0.5h}]

plot(p_ez, p_adv, layout=l) # dpi=900)
# begin
# save_loc_double = "plots\\NaI(Ti)_Double_ResponseTime_Fluoresence_Decay2.pdf"
# savefig(save_loc_double)
# end


## PMT fit

##
# begin
#     plot(pmtt, pmth, color=:blue, label="Collected data");
#     idx = 335
#     scatter!([pmtt[idx]], [pmth[idx]])
# end
# N_pmt_2(t, p) = @. p[3] * exp(-t/p[1]) + p[4] * exp(-t/p[2]) + p[5]/(p[6] - p[1]) * ( exp(-t/p[1]) - exp(-t/p[6]) )
# fit_pmt2 = curve_fit(N_pmt_2, ustrip.(t), ustrip.(h), p0pmt, lower=lbpmt)
# pmt_params2 = fit_pmt1.param

# plot!(pmtt, N_pmt_2(ustrip.(pmtt), ustrip.(pmt_params2)), label= "Fit to data, decay time: $(round((pmt_params2[1]),digits=2)*u"ns")",)

pmtt = df[!, :Time];
pmth = df[!, :Height];
N_pmt_2(t, p) = @. p[3] * exp(-t/p[1]) + p[4] * exp(-t/p[2]) + p[5]/(p[6] - p[1]) * ( exp(-t/p[1]) - exp(-t/p[6]) )

begin
    p0pmt = [230.0, 15.0, 120.0, 20.0, 1e-10, 5.0];
    lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-12, 1e-2];
    rise_start = 377
    fit_pmt2 = curve_fit(N_pmt_2, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), p0pmt, lower=lbpmt)
    fit_pmt2_start = curve_fit(N_pmt_2, ustrip.(pmtt[2:end]), ustrip.(pmth[2:end]), p0pmt, lower=lbpmt)
    midpoint = 227
    fit_pmt2_start_mid = curve_fit(N_pmt_2, ustrip.(pmtt[midpoint:end]), ustrip.(pmth[midpoint:end]), p0pmt, lower=lbpmt)
    pmt_params2 = fit_pmt2.param
    pmt_params2_start = fit_pmt2_start.param
    pmt_params2_start_mid = fit_pmt2_start_mid.param
    vars = ["τ_s", "τ_f", "A", "B", "GNeR", "τ"]
    @show vars
    @show round.(pmt_params2; digits=2)

    plot(pmtt, pmth, color=:blue, label="Collected data");
    plot!(pmtt[rise_start:end], N_pmt_2(ustrip.(pmtt[rise_start:end]), ustrip.(pmt_params2)),
            color=:red, label= "Fit to data, decay time: $(round((pmt_params2[1]), digits=2)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7,
            axisfontsize=4
    ) #
    plot!(pmtt[2:end], N_pmt_2(ustrip.(pmtt[2:end]), ustrip.(pmt_params2_start)),
            color=:green, label= "Fit to data (from start), decay time: $(round((pmt_params2_start[1]), digits=2)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7,
            axisfontsize=4,
            linstyle=:dot
    ) #
    plot!(pmtt[midpoint:end], N_pmt_2(ustrip.(pmtt[midpoint:end]), ustrip.(pmt_params2_start_mid)),
            label= "Fit to data (from $(midpoint)ns), decay time: $(round((pmt_params2_start_mid[1]), digits=2)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7,
            axisfontsize=4, color=:darkblue,
            linstyle=:dot
    ) #
    scatter!([pmtt[rise_start]], [pmth[rise_start]], color=:darkred, label="Plot from which i start the red fit")
    scatter!([pmtt[midpoint]], [pmth[midpoint]], label="Plot from which i start the purple fit", color=:orange)
end
savefig("plots\\help_me.png")
