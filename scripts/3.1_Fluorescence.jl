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

##
pmtt = df[!, :Time];
pmth = df[!, :Height];
N_pmt(t, p) = @. p[1] * exp(-t/p[2]) + p[3]*exp(-t/p[4]) - p[5] * exp(-t/p[6])

p_pmt = begin
    p0pmt = [120.0, 230.0, 80.0, 5.0, 50, 5.0];
    # lbpmt = zeros(length(p0pmt))
    lbpmt = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2];

    rise_start = 325
    fit_pmt = curve_fit(N_pmt, ustrip.(pmtt[rise_start:end]), ustrip.(pmth[rise_start:end]), p0pmt, lower=lbpmt)
    pmt_params = fit_pmt.param
    # pmt_errors = confidence_interval(fit_pmt)
    # @show round.(pmt_params; digits=2)

    pmt_fit = plot(pmtt, pmth, color=:blue, label="Collected data");
    plot!(pmtt[rise_start:end], N_pmt(ustrip.(pmtt[rise_start:end]), ustrip.(pmt_params)),
            color=:red, label= "Fit to data",#, decay time: $(round((pmt_params[1]), digits=2)*u"ns")",
            xlabel="Time series",
            ylabel="Height (mV)",
            legendfontsize=7,
            axisfontsize=4
    ) #
    latstring = let p = round.(pmt_params; digits=2)
        A = p[1]
        τs = p[2]
        B = p[3]
        τf = p[4]
        C = p[5]
        τ = p[6]
        (L"""A = %$(A), \tau_s = %$(τs)""", L"""B = %$(B), \tau_f = %$(τf)""", L"""C=%$(C), \tau=%$(τ)""")
    end

    scatter!([pmtt[rise_start]], [pmth[rise_start]], color=:darkred, label="Plot from which i start the red fit")
    annotate!(750, 85, text(L"""A e^{-\frac{t}{\tau_s}} + B e^{-\frac{t}{\tau_f}} - C e^{-\frac{t}{\tau}}""", :red, 12));
    annotate!(750, 72, text(latstring[1], 10) );
    annotate!(750, 62, text(latstring[2], 10) );
    annotate!(750, 52, text(latstring[3], 10) );
end
# savefig("plots\\help_me.png")

l = @layout [a{0.5h};
             b{0.5h}]
plot(p_ez, p_pmt, layout=l)
