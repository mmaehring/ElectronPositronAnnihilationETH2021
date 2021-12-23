## Marcus Mähring 2021

#=
This program fits the gaussian peaks and compton edges which appear in the MCA
recorded energy spectrum. The idea is to
1) Find the corresponding energies for the edges and the peaks
2) Perform a linear fit to determine a conversion from ADC to keV
3) Determine if this fit is indeed linear -> χ² test
4) Find energy of low energy peak


The fits to the peaks will be done using a gaussian fit with different background.
One should note that the amount of sources contributing to the structure of
the spectrum decreases with higher energies. Thus, we start with fitting the first peak with
constant noise from the detector.
Then for the compton edge we also account for linear noise
For the photopeak at 511 keV we determine the fit using quadratic noise
The same noise is used at the compton edge, because I'm too stupid to come up with
a better alternative.
=#
## Loading packages
begin
    using DrWatson, Logging
    # Plotting
    using StatsPlots, UnitfulRecipes, Measures, LaTeXStrings
    # Dealing with data
    using Unitful, Measurements, DataFrames, CSV, Statistics, StatsBase
    # Fitting
    using LsqFit, Distributions
    # additional utilities
    using SpecialFunctions, Interpolations, Combinatorics
end

## Helper functions
function χ²_calc(fit_function, x, y, σ, dependent_vars = 3)
        χ² = sum((fit_function.(x) .- y).^2 ./ σ.^2)
        dof = length(x) - dependent_vars - 1
        χ²dof = χ²/dof
        return dof, χ², χ²dof
end

function χ²_test(χ²_value, degrees_of_freedom; α=0.05)
    #=
    H₀: There is no reason to believe that the model describes the data than any other
    If χ² < cutoff -> reject null hypothesis at confidence α
    =#
    # @info degrees_of_freedom
    χ²_distr = Distributions.Chisq(degrees_of_freedom)
    χ²_cutoff = Distributions.quantile(χ²_distr, 1 - α)
    p_val_χ² = (1 - Distributions.cdf(χ²_distr, χ²_value))
    best_possible = Distributions.cdf(χ²_distr, χ²_value)
    return χ²_cutoff, χ²_value, χ²_value >= χ²_cutoff, p_val_χ², best_possible, χ²_value <= χ²_cutoff
end

function compute_runs(data, model_data)
    arr = 0 .> (model_data .- data)
    current = arr[1]
    rtmp = 1
    i = 1
    while i != length(arr)
        i = i+1
        if arr[i] != current
            rtmp += 1
            current = arr[i]
        end
    end
    rtmp
end

# Less than 10-15 points -> otherwise gaussian.
function run_test(data, model_data)
    N = length(data)
    Na = sum(0 .> (model_data .- data))
    Nb = sum(0 .< (model_data .- data))

    @debug 0 .< (model_data .- data)
    @debug 0 .> (model_data .- data)

    r_data = compute_runs(data, model_data)

    @debug r_data

    # Nb = N - Na
    if N != Nb + Na
        @warn "Abandon printf debugging, all ye who enter here!: Na + Nb != N, some points are on the line + ϵ"
    end

    r_exp = 1 + 2 * Na * Nb / N
    V(r) = 2*Na*Nb*(2 * Na * Nb - N) / (N^2 * (N - 1))
    σ(r) = sqrt(V(r))

    return r_data, r_exp, V(r_exp), σ(r_exp), N, Na
end

function prob_r_runs(data, model_data)
    N = length(data)
    Na = sum(0 .> (model_data .- data))
    Nb = sum(0 .< (model_data .- data))

    @debug 0 .> (model_data .- data)
    @debug 0 .< (model_data .- data)

    r_data = compute_runs(data, model_data)
    if r_data % 2 == 0
         tmp = binomial(Na-1, (0.5*r_data - 1) |> Int64) * binomial(Nb-1, (0.5*r_data - 1) |> Int64) / binomial(N, Na)
         @debug tmp
        return 2*tmp;
    else
        @debug (r_data-1)/2 == (r_data-1)/2 |> Int64

        tmp = binomial(Na-1, (r_data - 3)/2 |> Int64) * binomial(Nb-1, (r_data-1)/2 |> Int64) +
              binomial(Na-1, (r_data-1)/2 |> Int64) * binomial(Nb-1, (r_data-3)/2 |> Int64)
        tmp /= binomial(N, Na)
        @debug tmp
        return tmp;
    end
end

function poisson_coinfidence_band(ADC, Counts; α = 0.5, σC = 2)
    CountError = sqrt.(Counts)
    l = (Counts .- σC*CountError)
    u = (Counts .+ σC*CountError)
    c_confidence=:red
    plot!(ADC, Counts, fillrange = u, linealpha=0,
          fillalpha = α, c = c_confidence, label = "Confidence band ($(σC)σ)", st=:step, linewidth=0
    )
    plot!(ADC, Counts, fillrange = l, st=:step,
          fillalpha = α, c = c_confidence, label = :none, linewidth=0, linealpha=0,
    )
end


## Defining the models

gaussian_f(x,p) = @.p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) )
gaussian_quad_noise_f(x,p) = @. p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) ) + p[4] + p[5]*x + p[6]*x^2;
gaussian_lin_noise_f(x,p) = @. p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) ) + p[4] + p[5]*x;
gaussian_const_noise_f(x,p) = @. p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) ) + p[4];

errorfunc_normal_gauss_const_noise_f(x,p) = @. p[1]*erf(x)+ p[2];
errorfunc_arb_gauss_const_noise_f(x,p) = @. p[1] * erf( (x-p[2]) / (sqrt(2)*p[3]) ) + p[4];
errorfunc_arb_gauss_linear_noise_f(x,p) = @. p[1] * erf( (x-p[2]) / (sqrt(2)*p[3]) ) + p[4] + p[5]*x;
errorfunc_arb_gauss_double_linear_noise_f(x,p) = @. p[1] * erf( (x-p[2]) / (sqrt(2)*p[3]) ) + p[4] + p[5]*x + p[7] + p[8]*x;
errorfunc_arb_gauss_quad_noise_f(x,p) = @. p[1] * erf( (x-p[2]) / (sqrt(2)*p[3]) ) + p[4] + p[5]*x + p[6]*x^2 # didn't work

## Loading the data

path = "data\\3.3\\3.3.2_G\\3.3.2.txt"; # data is stored in 3.3.2 -> this was the acquisition
spectrum_data = CSV.read(path, DataFrame);
spectrum_data[!, :Error] = sqrt.(spectrum_data[!, :Count])


## Initial histogram plot of the data set

begin
    p_uncal_spectrum = plot(spectrum_data[!, :Channel], spectrum_data[!, :Count],
        seriestype=:step,
        xlims=(0,1000), xticks=0:100:1000, xlabel="Channel",
        title="MCA spectrum - Uncalibrated", ylabel="Detected counts",
        size=(700,500), label="Counts", dpi=700, titlefontsize=10, xguidefontsize=9, yguidefontsize=9,
        margin = 5mm
    )
    # poisson_coinfidence_band(spectrum_data[!, :Channel], spectrum_data[!, :Count], α=0.5)
    # savefig("plots\\MCA_Spectrum_NoCoinc_RAW_Errors.pdf")
    fontsize=7
    colors = distinguishable_colors(12)[5:end]
    # annotate!(80, 2.07e4, text("Electron capture", colors[1] , fontsize))
    annotate!(215, 5.1e3, text("511 keV Compton edge", colors[2] , fontsize, rotation = -70.0))
    # annotate!(120, 1.1e4, text("Backscattering",  colors[3], fontsize))
    annotate!(310, 2.175e4, text("511 keV peak",  colors[4], fontsize))
    annotate!(600, 1.7e3, text("1275 keV Compton edge",  colors[6], fontsize))
    annotate!(780, 2.5e3, text("1275 keV peak",  colors[8], fontsize))
    # savefig("plots\\3.3.2_MCA_Spectrum_NoCoinc_Annotated_Mod_LeftPeaksUnmarked.pdf")
end
## 1275 peak   # Frequentist  -> there is very little noise here -> no background fit

p0 = [2250.0, 780.0, 20.0]
p0_const = [2250.0, 780.0, 20.0, 100.0]

Peak1275_data = spectrum_data[725:825, :];
Peak1275_x_data = Peak1275_data[!, :Channel]
Peak1275_y_data = Peak1275_data[!, :Count]
Peak1275_σ_data = Peak1275_data[!, :Error]

Peak1275_fit = curve_fit(gaussian_f, Peak1275_x_data, Peak1275_y_data, 1 ./Peak1275_σ_data.^2 ,p0)
Peak1275_const_fit = curve_fit(gaussian_const_noise_f, Peak1275_x_data, Peak1275_y_data, 1 ./Peak1275_σ_data.^2, p0_const)
Peak1275_f(x) = gaussian_f(x , Peak1275_fit.param)
Peak1275_const_f(x) = gaussian_const_noise_f(x, Peak1275_const_fit.param)
Peak1275_no_noise_chisqndof = χ²_calc(Peak1275_f, Peak1275_x_data, Peak1275_y_data, Peak1275_y_data, length(p0))[end]
Peak1275_const_noise_chisqndof = χ²_calc(Peak1275_const_f, Peak1275_x_data, Peak1275_y_data, Peak1275_y_data, length(p0))[end]

values_1275 = Peak1275_fit.param
sigma1275 = stderror(Peak1275_fit) # cov_mat1275 = estimate_covar(Peak1275_fit)

values_const_1275 = Peak1275_const_fit.param;
sigma_const_1275 = stderror(Peak1275_const_fit);

begin
    p_1275= plot(spectrum_data[!, :Channel][720:825], spectrum_data[!, :Count][720:825],
        size=(800,500), seriestype=:step, title="1275 keV peak", xlabel="ADC", ylabel="Counts",
        label="Recorded data", margin=3mm, xlims=(710,827),
        legend=:topleft, titlefontsize=12
    )
    poisson_coinfidence_band(spectrum_data[!, :Channel][720:825], spectrum_data[!, :Count][720:825], α=0.15)
    plot!(720:825, Peak1275_f, label="Fit w.o. noise, χ²/ndof = $(round(Peak1275_no_noise_chisqndof, digits=3))", color=:darkblue)
    plot!(720:825, Peak1275_const_f, label="Fit w. constant noise, χ²/ndof = $(round(Peak1275_const_noise_chisqndof, digits=3))", color=:darkred)
    annotate!(812, 1900, text("μ = $(round(values_1275[2], digits=2)) ± $(round.(sigma1275[2], digits=2))", 10, :darkblue))
    annotate!(812, 1700, text("σ = $(round(values_1275[3], digits=2)) ± $(round.(sigma1275[3], digits=2))", 10, :darkblue))
    annotate!(812, 1500, text("μ = $(round(values_const_1275[2], digits=2)) ± $(round.(sigma_const_1275[2], digits=2))", 10, :darkred))
    annotate!(812, 1300, text("σ = $(round(values_const_1275[3], digits=2)) ± $(round.(sigma_const_1275[3], digits=2))", 10, :darkred))
    # savefig("plots\\3.3.3_1275_keV_Peak_w_sigma.pdf")
end

## 1275 edge
CE1275_ADC = spectrum_data[!, :Channel][550:725]
CE1275_CNT = spectrum_data[!, :Count][550:725]
CE1275_CNTσ = sqrt.(spectrum_data[!, :Count][550:725])
# corr_indx = first(findall(x->x==630, CE1275_ADC)) # ADC channel 625 in the new array -> about here the edge begins
corr_indx = first(findall(x->x==625, CE1275_ADC)) #
# corr_indx = first(findall(x->x==610, CE1275_ADC)) #

begin
    p01275CE = [-1000.0, 625.0, 60.0, 1100.0]
    CE1275_fit = curve_fit(errorfunc_arb_gauss_const_noise_f, (CE1275_ADC[corr_indx:end]), (CE1275_CNT[corr_indx:end]), 1 ./ (CE1275_CNTσ[corr_indx:end]).^2 , p01275CE)
    CE1275_f(x) = errorfunc_arb_gauss_const_noise_f(x , CE1275_fit.param);
    CE1275Coeffs = CE1275_fit.param .± stderror(CE1275_fit);
    CE1275CoeffsRounded = round.(CE1275_fit.param .± stderror(CE1275_fit), digits=2);
    CE1275Constχ = χ²_calc(CE1275_f, CE1275_ADC[corr_indx:end], CE1275_CNT[corr_indx:end], CE1275_CNTσ[corr_indx:end], length(p01275CE))[end]
end

begin
    p01275CE_linear = [-1200.0, 625.0, 60, 1100.0, 5.0]
    CE1275_linear_fit = curve_fit(errorfunc_arb_gauss_linear_noise_f, (CE1275_ADC[corr_indx:end]), (CE1275_CNT[corr_indx:end]), 1 ./ (CE1275_CNTσ[corr_indx:end]).^2 ,  p01275CE_linear)
    CE1275_linear_f(x) = errorfunc_arb_gauss_linear_noise_f(x, CE1275_linear_fit.param)
    CE1275LinCoeffs = CE1275_linear_fit.param .± stderror(CE1275_linear_fit);
    CE1275LinCoeffsRounded = round.(CE1275_linear_fit.param .± stderror(CE1275_linear_fit), digits=2);
    CE1275Linχ = χ²_calc(CE1275_linear_f, CE1275_ADC[corr_indx:end], CE1275_CNT[corr_indx:end], CE1275_CNTσ[corr_indx:end], length(p01275CE))[end]
end

begin # Here CE stands for compton edge
    p_compton_edge_1275 = plot(CE1275_ADC, CE1275_CNT, color=1,
        seriestype=:step, size=(750,400), label="Recorded data", margin=3mm,
        titlefontsize=10, dpi=900, xlabel="ADC", ylabel="Counts", xguidefontsize=9, yguidefontsize=9,#ylims=(0, 2500),
        title="Compton edge of 1275keV", legend=:bottomleft, legendfontsize=7
    )
    poisson_coinfidence_band(CE1275_ADC, CE1275_CNT, α=.15)
    # Set fontsize
    fontsize=8
    l_w = 0.8 # line width of the plotted lines

    # Model we're fitting to
    annotate!(685, 1700, text(L"A \cdot erf{ \left( \frac{x - \mu}{\sqrt{2} \sigma} \right)} + b + k \cdot x", :black, 11))

    # Plotting constant solution
    c1 = :purple
    plot!(CE1275_ADC, CE1275_f, label="Fit, χ²/ndof = $(round(CE1275Constχ, digits=2))", color=:purple, linewidth=l_w)
    annotate!(685, 1500, text("A = $(CE1275CoeffsRounded[1])", c1, fontsize))
    annotate!(685, 1400, text("μ = $(CE1275CoeffsRounded[2])", c1, fontsize))
    annotate!(685, 1300, text("σ = $(CE1275CoeffsRounded[3])", c1, fontsize))
    annotate!(685, 1200, text("b = $(CE1275CoeffsRounded[4])", c1, fontsize))
    annotate!(685, 1100, text("k = $(0)", c1, fontsize))
    vline!([Measurements.value(CE1275CoeffsRounded[2])], label="Compton edge w.o. scaling", color=c1, linestyle=:dash)

    # Plotting linear scaling solution
    c2 = :green
    plot!(CE1275_ADC, CE1275_linear_f, label="Fit with scaling, χ²/ndof = $(round(CE1275Linχ, digits=2))", color=:green, linewidth=l_w)
    annotate!(685, 950, text("A = $(CE1275LinCoeffsRounded[1])", c2,  fontsize))
    annotate!(685, 850, text("μ = $(CE1275LinCoeffsRounded[2])", c2, fontsize))
    annotate!(685, 750, text("σ = $(CE1275LinCoeffsRounded[3])", c2, fontsize))
    annotate!(685, 650, text("b = $(CE1275LinCoeffsRounded[4])", c2, fontsize))
    annotate!(685, 550, text("k = $(CE1275LinCoeffsRounded[5])", c2, fontsize))
    vline!([Measurements.value(CE1275LinCoeffsRounded[2])], label="Compton edge w. scaling", color=c2, linestyle=:dash)
    # savefig("plots\\3.3.3_1275CE.pdf")
end


## 1275 Compton continuum fitting
# quadMod(x,p) = @. p[3]*x^2 + p[2]*x + p[1]
# fitcomptoncontinuum = curve_fit(quadMod,
#                                spectrum_data[!, :Channel][350:600], spectrum_data[!, :Count][350:600], (spectrum_data[!, :Error][350:600]).^(-2) ,
#                                [20.0, 30.0, 1400.0, 2300.0]
# )
#
# p_compton_cont_1275 = plot(spectrum_data[!, :Channel][350:700], spectrum_data[!, :Count][350:700],
#     seriestype=:step,
#     title="Compton continuum and edge of 1275keV"
# )
# plot!(spectrum_data[!, :Channel][350:700], quadMod(spectrum_data[!, :Channel][350:700], fitcomptoncontinuum.param))


## 511 Peak
p0 = [22500.0, 325, 25.0]
p0_quad = [22500.0, 325, 25.0, 2000.0, 1, 0.1]
p0_lin = [22500.0, 325, 25.0, 2000.0, 0.5]
p0_const = [22500.0, 325, 25.0, 2000.0]

Peak511_data = spectrum_data[250:350, :];
Peak511_x_data = Peak511_data[!, :Channel];
Peak511_y_data = Peak511_data[!, :Count];
Peak511_σ_data = Peak511_data[!, :Error];

function gaussian_fit_to_model(gaussian_mod, xdata, ydata, σdata, p0_val)
    gaussian__fit = curve_fit(gaussian_mod, xdata, ydata, σdata.^-2, p0_val);
    Peak__f(x) = gaussian_f(x , gaussian__fit.param);
    try
        Peak__P = gaussian__fit.param .± stderror(gaussian__fit)
    catch e
        @error "Could not calculate error due to: $(e)"
        return -1
    end
    Peak__P = gaussian__fit.param .± stderror(gaussian__fit)
    Peak__quality  = χ²_calc(Peak__f, xdata, ydata,
                              σdata, length(gaussian__fit.param))
    return gaussian__fit, Peak__f, Peak__P, Peak__quality
end

# NN -> No noise
G511NNFit, G511NN_f, G511NN_peak, QNN = gaussian_fit_to_model(gaussian_f, Peak511_x_data, Peak511_y_data,
                                      Peak511_σ_data, p0)

begin
    gaussian_511_fit = curve_fit(gaussian_f, Peak511_x_data, Peak511_y_data,  Peak511_σ_data.^-2, p0);
    Peak511_P = gaussian_511_fit.param .± stderror(gaussian_511_fit)
    Peak511_f(x) = gaussian_f(x , gaussian_511_fit.param);
    Peak511_quality_R  = χ²_calc(Peak511_f, Peak511_x_data, Peak511_y_data,
                                Peak511_σ_data, length(gaussian_511_fit.param))
end

begin
    gaussian_quad_511_fit = curve_fit(gaussian_quad_noise_f, Peak511_x_data, Peak511_y_data, Peak511_σ_data.^-2, p0_quad);
    Peak511_quad_P = gaussian_quad_511_fit.param .± stderror(gaussian_quad_511_fit)
    Peak511_quad_f(x) = gaussian_quad_noise_f(x , gaussian_quad_511_fit.param);
    Peak511_quad_quality_R  = χ²_calc(Peak511_quad_f, Peak511_x_data, Peak511_y_data,
                                    Peak511_σ_data, length(gaussian_quad_511_fit.param))
end

# begin
#     heatmap(log.(abs.(estimate_covar(gaussian_quad_511_fit)[:,:])), c=:rainbow);
#     for (ind, val) in enumerate(["A", "μ", "σ", "c", "b", "a"])
#         annotate!(ind, ind, text(val, :black, 12));
#     end
#     plot!()
# end

begin
    gaussian_lin_511_fit = curve_fit(gaussian_lin_noise_f, Peak511_x_data, Peak511_y_data, Peak511_σ_data.^-2, p0_lin);
    Peak511_lin_P = gaussian_lin_511_fit.param .± stderror(gaussian_lin_511_fit)
    Peak511_lin_f(x) = gaussian_lin_noise_f(x, gaussian_lin_511_fit.param);
    Peak511_lin_quality_R  = χ²_calc(Peak511_lin_f, Peak511_x_data, Peak511_y_data,
                                    Peak511_σ_data, length(gaussian_lin_511_fit.param))
end

begin
    gaussian_const_511_fit = curve_fit(gaussian_const_noise_f, Peak511_x_data, Peak511_y_data, Peak511_σ_data.^-2, p0_const);
    Peak511_const_P = gaussian_const_511_fit.param .± stderror(gaussian_const_511_fit)
    Peak511_const_f(x) = gaussian_const_noise_f(x, gaussian_const_511_fit.param);
    Peak511_const_quality_R  = χ²_calc(Peak511_const_f, Peak511_x_data, Peak511_y_data,
                                       Peak511_σ_data, length(gaussian_const_511_fit.param))
end


begin
    p_511= plot(spectrum_data[!, :Channel], spectrum_data[!, :Count],
        seriestype=:step, label="Recorded data",
        xlims=(250, 350), size=(800,500), xlabel="ADC", ylabel="Count",
        title="511keV peak", legend=:topleft, margin=3mm
    )
    poisson_coinfidence_band((spectrum_data[!, :Channel][250:350]), (spectrum_data[!, :Count][250:350]);
        α=0.15, σC = 2
    )
    plot!(250:350, Peak511_quad_f, label="Quadractic: χ²/ndof = $(round(Peak511_quad_quality_R[end], digits=3))", c=:darkblue)
    plot!(250:350, Peak511_lin_f, label="Linear: χ²/ndof = $(round(Peak511_lin_quality_R[end], digits=3))", c=:darkred)
    plot!(250:350, Peak511_const_f, label="Constant: χ²/ndof = $(round(Peak511_const_quality_R[end], digits=3))", c=:cyan)
    plot!(250:350, Peak511_f, label="No removal: χ²/ndof = $(round(Peak511_quality_R[end], digits=3))", c=:darkgreen)
    annotate!(337, 1.95e4, text("μ=$(round(Peak511_quad_P[2], digits=3))0", :darkblue, 10))
    annotate!(337, 1.80e4, text("σ=$(round(Peak511_quad_P[3], digits=3))", :darkblue, 10))
    annotate!(337, 1.65e4, text("μ=$(round(Peak511_lin_P[2], digits=3))", :darkred, 10))
    annotate!(337, 1.50e4, text("σ=$(round(Peak511_lin_P[3], digits=3))", :darkred, 10))
    annotate!(337, 1.35e4, text("μ=$(round(Peak511_const_P[2], digits=3))", :cyan, 10))
    annotate!(337, 1.20e4, text("σ=$(round(Peak511_const_P[3], digits=3))", :cyan, 10))
    annotate!(337, 1.05e4, text("μ=$(round(Peak511_const_P[2], digits=3))", :darkgreen, 10))
    annotate!(337, .9e4, text("σ=$(round(Peak511_const_P[3], digits=3))", :darkgreen, 10))
    # savefig("plots\\3.3.3_511_keV_Peak_w_sigma.pdf")
end



## 511 edge
p_compton_continuum_511 = plot(spectrum_data[!, :Channel], spectrum_data[!, :Count],
    seriestype=:step,
    xlims=(0, 250),
    xticks=0:50:250,
    title="Compton continuum and edge of 511keV"
)

CE511_ADC = spectrum_data[!, :Channel][125:250];
CE511_CNT = spectrum_data[!, :Count][125:250];
CE511_CNTσ = sqrt.(spectrum_data[!, :Count][125:250]);

# corr_indx = first(findall(x->x==150, CE511_ADC))
corr_indx = first(findall(x->x==150, CE511_ADC))

begin
    p0511CE= [-1000.0, 175.0, 60, 1100.0]
    CE511_fit = curve_fit(errorfunc_arb_gauss_const_noise_f,
        (CE511_ADC[corr_indx:end]), (CE511_CNT[corr_indx:end]),
        1 ./ (CE511_CNTσ[corr_indx:end]).^2 , p0511CE
    );
    CE511_f(x) = errorfunc_arb_gauss_const_noise_f(x , CE511_fit.param);
    CE511Coeffs = CE511_fit.param .± stderror(CE511_fit);
    CE511CoeffsRounded = round.(CE511Coeffs, digits=2);
    χ²ndof_CE511_const = χ²_calc(CE511_f, CE511_ADC[corr_indx:end], CE511_CNT[corr_indx:end],
                              CE511_CNTσ[corr_indx:end], length(p0511CE))
end

begin
    p0511_lin_CE= [-1000.0, 175.0, 60, 1100.0, 5.0];
    CE511_linear_fit = curve_fit(errorfunc_arb_gauss_linear_noise_f,
        (CE511_ADC[corr_indx:end]), (CE511_CNT[corr_indx:end]),
        1 ./ (CE511_CNTσ[corr_indx:end]).^2 ,  p0511_lin_CE
    );
    CE511_linear_f(x) = errorfunc_arb_gauss_linear_noise_f(x, CE511_linear_fit.param);
    CE511LinCoeffs = CE511_linear_fit.param .± stderror(CE511_linear_fit);
    CE511LinCoeffsRounded = round.(CE511_linear_fit.param .± stderror(CE511_linear_fit), digits=2);
    χ²ndof_CE511_linear = χ²_calc(CE511_linear_f, CE511_ADC[corr_indx:end], CE511_CNT[corr_indx:end],
                              CE511_CNTσ[corr_indx:end], length(p0511_lin_CE))
end

# begin
#     p0511_lin_CE= [-1000.0, 175.0, 60, 1100.0, 5.0];
#     CE511_linear_fit = curve_fit(errorfunc_arb_gauss_linear_noise_f,
#         (CE511_ADC[corr_indx:end]), (CE511_CNT[corr_indx:end]),
#         1 ./ (CE511_CNTσ[corr_indx:end]).^2 ,  p0511_lin_CE
#     );
#     CE511_linear_f(x) = errorfunc_arb_gauss_linear_noise_f(x, CE511_linear_fit.param);
#     CE511LinCoeffs = CE511_linear_fit.param .± stderror(CE511_linear_fit);
#     CE511LinCoeffsRounded = round.(CE511_linear_fit.param .± stderror(CE511_linear_fit), digits=2);
#     χ²ndof_CE511_linear = χ²_calc(CE511_linear_f, CE511_ADC[corr_indx:end], CE511_CNT[corr_indx:end],
#                               CE511_CNTσ[corr_indx:end], length(p0511_lin_CE))
# end

begin # Here CE stands for compton edge
    p_compton_edge_511 = plot(CE511_ADC, CE511_CNT, color=1,
        seriestype=:step, size=(750,400), label=:none,
        titlefontsize=10, dpi=900, ylims = (1500, 10000), xlabel="ADC", ylabel="Count",
        title="Compton edge of 511keV", legend=:bottomleft, legendfontsize=8, margin=3mm
    )

    # scatter!([CE511_ADC[72]], [CE511_CNT[72]])

    # error band
    # l = (CE511_CNT .- CE511_CNTσ)
    # u = (CE511_CNT .+ CE511_CNTσ)
    # c_confidence = :red
    # α = 0.2
    # plot!(CE511_ADC, CE511_CNT, st=:step, fillrange = u, fillalpha = α, c = c_confidence, label = "Confidence band", linewidth=0)
    # plot!(CE511_ADC, CE511_CNT, st=:step, fillrange = l, fillalpha = α, c = c_confidence, label = :none, linewidth=0)

    # Set fontsize
    fontsize=9
    l_w = 1.5 # line width of the plotted lines
    # Model we're fitting to
    annotate!(220, 8750, text(L"A \cdot erf{ \left( \frac{x - \mu}{\sqrt{2} \sigma} \right)} + b + k \cdot x", :black, 11))
    c_confidence=:red

    # Plotting constant solution
    c1 = :darkblue
    plot!(CE511_ADC, CE511_f, label="Fit, χ²/ndof = $(round(χ²ndof_CE511_const[end], digits=2))", color=:darkblue, linewidth=l_w)
    annotate!(215, 7700, text("A = $(CE511CoeffsRounded[1])", c1, fontsize))
    annotate!(215.25, 7300, text("μ = $(CE511CoeffsRounded[2])", c1, fontsize))
    annotate!(215, 6900, text("σ = $(CE511CoeffsRounded[3])", c1, fontsize))
    annotate!(215.25, 6500, text("b = $(CE511CoeffsRounded[4])", c1, fontsize))
    annotate!(215, 6100, text("k = $(0)", c1, fontsize))
    vline!([Measurements.value(CE511CoeffsRounded[2])], label="Identified Compton edge", color=c1, ls=:dash)

    # Plotting linear scaling solution
    c2 = :darkred
    plot!(CE511_ADC, CE511_linear_f, label="Fit with linear scaling, χ²/ndof = $(round(χ²ndof_CE511_linear[end], digits=2))", color=:darkred, linewidth=l_w)
    annotate!(215, 5400, text("A = $(CE511LinCoeffsRounded[1])", c2,  fontsize))
    annotate!(215, 5000, text("μ = $(CE511LinCoeffsRounded[2])", c2, fontsize))
    annotate!(215, 4600, text("σ = $(CE511LinCoeffsRounded[3])", c2, fontsize))
    annotate!(215, 4200, text("b = $(CE511LinCoeffsRounded[4])", c2, fontsize))
    annotate!(215, 3800, text("k = $(CE511LinCoeffsRounded[5])", c2, fontsize))
    vline!([Measurements.value(CE511LinCoeffsRounded[2])], label="Identified Compton edge with scaling", color=c2, ls=:dash)

    poisson_coinfidence_band(CE511_ADC, CE511_CNT, α=0.15, σC = 2)

    # savefig("plots\\3.3.3_511CE.pdf")
end

# c3 = :orange
# plot!(CE511_ADC, errorfunc_arb_gauss_linear_noise_f(CE511_ADC, [-2500.0, 185.0, 16.0, 7850, -15.25, 0.6]), label="Fit with quadratic scaling ", color=c3, linewidth=l_w)
# plot!(CE511_ADC, errorfunc_arb_gauss_linear_noise_f(CE511_ADC, CE511_quad_fit.param), label="Fit with qudratic scaling ", color=:red, linewidth=l_w)
# # annotate!(215, 4600, text("A = $(CE511QuadCoeffsRounded[1])", c3,  fontsize))
# # annotate!(215.25, 4200, text("μ = $(CE511QuadCoeffsRounded[2])", c3, fontsize))
# # annotate!(215.05, 3800, text("σ = $(CE511QuadCoeffsRounded[3])", c3, fontsize))
# # annotate!(215.5, 3400, text("b = $(CE511QuadCoeffsRounded[4])", c3, fontsize))
# # annotate!(215, 3000, text("k = $(CE511QuadCoeffsRounded[5])", c3, fontsize))
# # vline!([Measurements.value(CE511_quad_fit.param[2])], label="Identified Compton edge with scaling", color=c3, ls=:dash)

## Energy calibration
best_position_CE511 = Measurements.value(CE511LinCoeffs[2]) ± 2*Measurements.uncertainty(CE511LinCoeffs[2])
best_position_PP511 = Measurements.value(Peak511_quad_P[2]) ± 2*Measurements.uncertainty(Peak511_quad_P[2])
best_position_CE1275 = Measurements.value(CE1275LinCoeffs[2]) ± 2*Measurements.uncertainty(CE1275LinCoeffs[2])
best_position_PP1275 = Peak1275_const_fit.param[2] ± 2*stderror(Peak1275_const_fit)[2]
energies = [340.666667, 511, 850, 1275] # keV
ADC_vals = [best_position_CE511, best_position_PP511, best_position_CE1275, best_position_PP1275]

# linear_model(x,p) = @. p[1] + x*p[2]
# begin
#     initial = [0.0, 120.0]
#     w = (Measurements.uncertainty.(ADC_vals)).^(-2)
#     calibration_fit = curve_fit(linear_model,
#                     Measurements.value.(energies), Measurements.value.(ADC_vals),
#                     w,
#                     initial
#     )
#     coeffs = calibration_fit.param .± stderror(calibration_fit)
#     calibration_f(x) = linear_model(x, calibration_fit.param)
#     scatter(energies, ADC_vals, xlabel="Energies (keV)", ylabel="ADC",
#         title="Energy ADC conversion", titlefontsize=10,
#         legend=:topleft, label="Data ± 2σ"
#     )
#     @info coeffs
#     plot!(200:50:1325, calibration_f, label="Fit")
#     annotate!(850, 225,
#               text("ADC = $(round(coeffs[2],digits=4)) ⋅ Energy - $(round(-1*coeffs[1], digits=2))", 9)
#     )
#     # @info w
#     # @info w .* (ADC_vals - calibration_f.(energies)).^2
#     # savefig("plots\\3.3.3_LinearFit.pdf")
# end

## New errors
linear_model(x,p) = @. p[1] + x*p[2]
begin
    best_position_CE511_ALT_ERR = Measurements.value(CE511LinCoeffs[2]) ± 15
    best_position_PP511_ALT_ERR = Measurements.value(Peak511_quad_P[2]) ± 2
    best_position_CE1275_ALT_ERR = Measurements.value(CE1275LinCoeffs[2]) ± 20
    best_position_PP1275_ALT_ERR = Peak1275_const_fit.param[2] ± 2
    ADC_vals = [best_position_CE511_ALT_ERR, best_position_PP511_ALT_ERR, best_position_CE1275_ALT_ERR, best_position_PP1275_ALT_ERR]
end

Calibfit = begin
    initial = [0.0, 120.0]
    w = (Measurements.uncertainty.(ADC_vals)).^(-2)
    calibration_fit = curve_fit(linear_model,
                    Measurements.value.(energies), Measurements.value.(ADC_vals),
                    w,
                    initial
    )
    coeffs = calibration_fit.param .± stderror(calibration_fit)
    calibration_f(x) = linear_model(x, calibration_fit.param)
    scatter(energies, ADC_vals, ylabel="ADC",
        title="Energy ADC conversion", titlefontsize=10,
        legend=:topleft, label="Data ± σ", dpi=750, xguidefontsize = 9, yguidefontsize = 9,
    )
    # @info ADC_Vals
    @info χ²_calc(calibration_f, energies,
        Measurements.value.(ADC_vals), Measurements.uncertainty.(ADC_vals),
        2
    )
    plot!(200:50:1325, calibration_f, label="Fit, χ² = $( round( χ²_calc(calibration_f, energies, Measurements.value.(ADC_vals), Measurements.uncertainty.(ADC_vals), 2)[2], digits=2))")
    annotate!(850, 225,
              text("ADC = $(round(coeffs[2],digits=4)) ⋅ Energy - $(round(-1*coeffs[1], digits=2))", 9)
    )
    # savefig("plots\\3.3.1_energy_new_errors.pdf")
end

begin
    residual_plot = plot(energies, (ADC_vals .- calibration_f.(energies)),
        seriestype=:scatter,
        label = "Residuals",
        xticks=350:100:1300,
        xlabel="Energies (keV)",
        ylabel = "Fit residuals", xguidefontsize = 9, yguidefontsize = 9,
        legend = :topleft,
        legendfontsize=7
    )
    plot!((x->0), [350, 1300], linestyle = :dash, label=:none, linecolor=:red)
end

begin
    l = @layout [a{0.6h};
    b{0.4h}]
    plot(Calibfit, residual_plot, layout = l)
    # savefig("plots\\3.3.3_energy_fit_w_residuals.pdf")
end

pull_plot = begin
    pulls = (ADC_vals .- calibration_f.(energies)) ./ Measurements.uncertainty.(ADC_vals)
    pull_plot = plot(energies, pulls,
        seriestype = :scatter,
        label = "Pulls", dpi=750,
        xlabel = "Energy [keV]", xguidefontsize = 9, yguidefontsize = 9,
        ylabel = "Pulls",
        legend = :top,
    )
    actual_mean_array = zeros(length(pulls)) .+ mean(pulls)
    # plot!(1:1:10, Measurements.value.(actual_mean_array), label = "Actual Average")
    plot!(energies, Measurements.value.(actual_mean_array), ribbon = zeros(4) .+ Measurements.uncertainty.(mean(pulls)) ,
        fillalpha = 0.25, c = 1, lw = 2,
        label = "Pull Average w. Error",
        linecolor = :red, legend=:topright,
        fillcolor = :red,
        legendfontsize = 6
    )
    plot!(energies, zeros(4), ribbon = ones(4) , fillalpha = 0.25, c = 1, lw = 2,
          label = "1σ around 0"
    )
end

hist_plt = begin
    histogram(Measurements.value.(pulls),
        weights = ones(length(pulls)) / length(pulls),
        xlabel = "Pulls",
        ylabel = "PDF",
        label = "Plot histogram"
    )

    fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(pulls))
    plot!(fitted_distribution, label="Fitted gaussian", legend=:topleft)
end
begin
    l_pull = @layout[a{0.75h}; b]
    plot(pull_plot, hist_plt, layout = l_pull)
    # savefig("plots\\3.3.3_energy_fit_pull_plot.pdf")
end
## fit quality

chisq_all_vals = χ²_calc(calibration_f, energies,
        Measurements.value.(ADC_vals), Measurements.uncertainty.(ADC_vals),
        2
)

χ²_test(chisq_all_vals[2], 1; α=0.05)

chisq_good_vals = χ²_calc(calibration_f, energies[BitArray([1,1,0,1])], Measurements.value.(ADC_vals[BitArray([1,1,0,1])]),
        Measurements.uncertainty.(ADC_vals[BitArray([1,1,0,1])]),
        2
)

χ²_test(chisq_good_vals[2], 1; α=0.05)
χ²_test(chisq_good_vals[2], 1)


run_test(ADC_vals, calibration_f(energies))
prob_r_runs(ADC_vals, calibration_f(energies))
# r_data, r_exp, V(r_exp), σ(r_exp), N, Na

## interlude on energy resolution
energy_to_adc(x) = calibration_f(x)
adc_to_energy(x) = (x - calibration_fit.param[1]) / calibration_fit.param[2]

FWHM511 = adc_to_energy.(Peak511_quad_P)[BitArray([0,1,1,0,0,0])]
RFWHM511 = 2.35 * FWHM511[2] / FWHM511[1]
RFWHM511_perc = 2.35 * FWHM511[2] / FWHM511[1] * 100

FWHM1275 = adc_to_energy.((Peak1275_const_fit.param .± stderror(Peak1275_const_fit))[BitArray([0,1,1,0])])
RFWHM1275 = 2.35 * FWHM1275[2] / FWHM1275[1]
RFWHM1275_perc = 2.35 * FWHM1275[2] / FWHM1275[1] * 100

## Fitting the peak at 0
begin
    PU_ADC = spectrum_data[!, :Channel][15:50];
    PU_CNT = spectrum_data[!, :Count][15:50];
    PU_CNTσ = sqrt.(spectrum_data[!, :Count][15:50]);
end

fitL, peak_f, pL, qL = gaussian_fit_to_model(gaussian_lin_noise_f, PU_ADC, PU_CNT, PU_CNTσ, [2.5e4, 27, 6, 1e4, 50])
fitQ, peak_Q, pQ, qQ = gaussian_fit_to_model(gaussian_quad_noise_f, PU_ADC, PU_CNT, PU_CNTσ, [2.5e4, 27, 6, 1e4, 50, 1.5])
fitC, peak_C, pC, qC = gaussian_fit_to_model(gaussian_const_noise_f, PU_ADC, PU_CNT, PU_CNTσ, [2.5e4, 27, 6, 1e4])
# gaussian_fit_to_model()

begin
    p_1275= plot(PU_ADC, PU_CNT,
        size=(800,500), seriestype=:step, title="Peak at low energy", xlabel="Energy (keV)", ylabel="Counts",
        label="Recorded data", margin=3mm, legend=:topright, titlefontsize=12
    )
    x_new = 14:0.1:45;
    poisson_coinfidence_band(PU_ADC, PU_CNT, α=0.15)
    plot!(x_new, gaussian_const_noise_f(x_new, Measurements.value.(pC)), label="Fit w. constant noise, χ²/ndof = $(round(qQ[end], digits=3))", color=:darkred)
    plot!(x_new, gaussian_lin_noise_f(x_new, Measurements.value.(pL)), label="Fit w. linear noise, χ²/ndof = $(round(qL[end], digits=3))", color=:darkblue)
    plot!(x_new, gaussian_quad_noise_f(x_new, Measurements.value.(pQ)), label="Fit w. quadratic noise, χ²/ndof = $(round(qQ[end], digits=3))", color=:darkgreen)
    annotate!(37.5, 16000, text("μ = $(round(adc_to_energy(pC[2]), digits=2))", 10, :darkred))
    annotate!(46, 16000, text("σ = $(round(adc_to_energy(pC[3]), digits=2))", 10, :darkred))
    annotate!(37.5, 15000, text("μ = $(round(adc_to_energy(pL[2]), digits=2))", 10, :darkblue))
    annotate!(46, 15000, text("σ = $(round(adc_to_energy(pL[3]), digits=2))", 10, :darkblue))
    annotate!(37.5, 14000, text("μ = $(round(adc_to_energy(pQ[2]), digits=2))", 10, :darkgreen))
    annotate!(46, 14000, text("σ = $(round(adc_to_energy(pQ[3]), digits=2))", 10, :darkgreen))
    plot!(xticks = (20:10:50, [string(round(adc_to_energy(i), digits=2)) for i in 20:10:52]))
    # savefig("plots\\3.3.3_low_keV_Peak.pdf")
end

# return gaussian__fit, Peak__f, Peak__P, Peak__quality
## Calibrated

begin
    p_compton_continuum_511 = plot(spectrum_data[!, :Channel], spectrum_data[!, :Count],
        seriestype=:step, label="Recorded data",
        xlims=(0, 250),
        xticks=(0:50:425, [string(round(adc_to_energy(i), digits=2)) for i in 0:50:425]),
        xlabel="Energy (keV)", xguidefontsize=9, yguidefontsize=9, titlefontsize=10,
        ylabel="Counts", lw=1.75, margin=3mm, dpi=750,
        title="Compton continuum and edge of 511keV"
    )
    poisson_coinfidence_band(spectrum_data[!, :Channel], spectrum_data[!, :Count], α=0.5, σC = 2)
    savefig("plots\\3.3.3calibrated_511_Compton_continuum.pdf")
end
