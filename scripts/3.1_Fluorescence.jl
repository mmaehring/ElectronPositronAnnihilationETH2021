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
path = "data\\3.1\\Dynode3.1\\800V\\800V_ALL0000_CurveData.txt"
begin
    df = CSV.read(path, DataFrame; select=[1,2]);
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
N_adv(t, p) = @. p[1] * exp(-t/p[2]) + p[3] * exp(-t/p[4]);

p0 = [500.0, 230.0, 500.0, 100.0];

fit2 = curve_fit(N_adv, ustrip.(t), ustrip.(h), p0);
adv_params = fit2.param
adv_errrors = stderror(fit2)

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

## plotta tillsammans

l = @layout [a{0.5h};
             b{0.5h}]

plot(p_ez, p_adv, layout=l) # dpi=900)

save_loc_double = "plots\\NaI(Ti)_Double_ResponseTime_Fluoresence_Decay2.pdf"
# savefig(save_loc_double)
