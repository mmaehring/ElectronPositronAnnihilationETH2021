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
    using DrWatson
    # Plotting
    using StatsPlots, UnitfulRecipes, Measures
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
        dof = length(x) - dependent_vars
        χ²dof = χ²/dof
        return dof, χ², χ²dof
end

function χ²_test(χ²_value, degrees_of_freedom; α=0.05)
    #=
    H₀: There is no reason to believe that the model describes the data than any other
    If χ² < cutoff -> reject null hypothesis at confidence α
    =#
    χ²_distr = Distributions.Chisq(degrees_of_freedom - 1)
    χ²_cutoff = Distributions.quantile(χ²_distr, 1 - α)
    return χ²_cutoff, χ²_value, χ²_value < χ²_cutoff
end

function ratio_test(data, model_data)
    

    return 0;
end

## Defining the models

gaussian_f(x,p) = @.p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) )
gaussian_quad_noise_f(x,p) = @. p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) ) + p[4] + p[5]*x + p[6]*x^2
gaussian_lin_noise_f(x,p) = @. p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) ) + p[4] + p[5]*x
gaussian_const_noise_f(x,p) = @. p[1] * exp( -(x - p[2])^2 / ( 2*p[3]^2) ) + p[4]

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
        size=(700,500), label="Counts", dpi=700, titlefontsize=11, xguidefontsize=10, yguidefontsize=10,
        margin = 5mm
    )
    # savefig("plots\\MCA_Spectrum_NoCoinc_RAW_Mod.pdf")

    #TODO -> fix the annotations, they are most likely wrong. Include electron capture and photopeaks = photoelectric effect Done???
    fontsize=7
    colors = distinguishable_colors(12)[5:end]
    # annotate!(80, 2.07e4, text("Electron capture", :green, fontsize))
    # annotate!(215, 5.1e3, text("511 keV Compton edge", :black, fontsize, rotation = -70.0))
    # annotate!(120, 1.1e4, text("Backscattering", :red, fontsize))
    # annotate!(310, 2.175e4, text("511 keV peak", :darkblue, fontsize))
    # annotate!(600, 1.7e3, text("1275 keV Compton edge", :purple, fontsize))
    # annotate!(780, 2.5e3, text("1275 keV peak", :orange, fontsize))
    annotate!(80, 2.07e4, text("Electron capture", colors[1] , fontsize))
    annotate!(215, 5.1e3, text("511 keV Compton edge", colors[2] , fontsize, rotation = -70.0))
    annotate!(120, 1.1e4, text("Backscattering",  colors[3], fontsize))
    annotate!(310, 2.175e4, text("511 keV peak",  colors[4], fontsize))
    annotate!(600, 1.7e3, text("1275 keV Compton edge",  colors[6], fontsize))
    annotate!(780, 2.5e3, text("1275 keV peak",  colors[8], fontsize))

    # savefig("plots\\MCA_Spectrum_NoCoinc_Annotated_Mod.pdf")
end
# p_uncal_log = plot(spectrum_data[!, :Channel], log.(10, spectrum_data[!, :Count]),
#     seriestype=:step,
#     xlims=(0,950),
#     xticks=0:100:900,
#     title="Entire spectrum, log",
#     yticks=0:0.5:5, label="Counts"
# )

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

values_1275 = Peak1275_fit.param
sigma1275 = stderror(Peak1275_fit); # cov_mat1275 = estimate_covar(Peak1275_fit);
# margin_of_error1275 = margin_error(Peak1275_fit, 0.05);
# confidence_inter1275 = confidence_interval(Peak1275_fit, 0.05);

values_const_1275 = Peak1275_const_fit.param
sigma_const_1275 = stderror(Peak1275_const_fit); # cov_mat1275 = estimate_covar(Peak1275_fit);
# margin_of_error1275 = margin_error(Peak1275_fit, 0.05);
# confidence_inter1275 = confidence_interval(Peak1275_fit, 0.05);

p_1275= plot(spectrum_data[!, :Channel][720:825], spectrum_data[!, :Count][720:825],
            seriestype=:step, title="1275 keV peak", label="Recorded data", legend=:topleft, legend
)
plot!(720:825, Peak1275_f, label="μ = $(round(values_1275[2], digits=2)) ± $(round.(sigma1275[2], digits=2))")
plot!(720:825, Peak1275_const_f, label="W noise, μ = $(round(values_const_1275[2], digits=2)) ± $(round.(sigma_const_1275[2], digits=2))")


## 1275 edge
p_compton_edge_1275 = plot(spectrum_data[!, :Channel][550:725], spectrum_data[!, :Count][550:725],
    seriestype=:step,
    xlims=(550, 725),
    title="Compton edge of 1275keV"
)


# ## 1275 Compton continuum fitting
# p_compton_cont_1275 = plot(spectrum_data[!, :Channel][350:700], spectrum_data[!, :Count][350:700],
#     seriestype=:step,
#     title="Compton continuum and edge of 1275keV"
# )
#

## 511 Peak
p0 = [22500.0, 325, 25.0]
p0_quad = [22500.0, 325, 25.0, 2000.0, 1, 0.1]
p0_lin = [22500.0, 325, 25.0, 2000.0, 0.5]
p0_const = [22500.0, 325, 25.0, 2000.0]

Peak511_data = spectrum_data[250:350, :];
Peak511_x_data = Peak511_data[!, :Channel];
Peak511_y_data = Peak511_data[!, :Count];
Peak511_σ_data = Peak511_data[!, :Error];

gaussian_511_fit = curve_fit(gaussian_f, Peak511_x_data, Peak511_y_data,  Peak511_σ_data.^-2, p0);
Peak511_P = gaussian_511_fit.param .± stderror(gaussian_511_fit)
Peak511_f(x) = gaussian_f(x , gaussian_511_fit.param);
Peak511_quality_R  = quality_of_fit_f(Peak511_f, Peak511_x_data, Peak511_y_data,
                                      Peak511_σ_data, length(gaussian_511_fit.param))

gaussian_quad_511_fit = curve_fit(gaussian_quad_noise_f, Peak511_x_data, Peak511_y_data, Peak511_σ_data.^-2, p0_quad);
Peak511_quad_P = gaussian_quad_511_fit.param .± stderror(gaussian_quad_511_fit)
Peak511_quad_f(x) = gaussian_quad_noise_f(x , gaussian_quad_511_fit.param);
Peak511_quad_quality_R  = quality_of_fit_f(Peak511_quad_f, Peak511_x_data, Peak511_y_data,
                                      Peak511_σ_data, length(gaussian_quad_511_fit.param))
begin
    heatmap(log.(abs.(estimate_covar(gaussian_quad_511_fit)[:,:])), c=:rainbow);
    for (ind, val) in enumerate(["A", "μ", "σ", "c", "b", "a"])
        annotate!(ind, ind, text(val, :black, 12));
    end
    plot!()
end

gaussian_lin_511_fit = curve_fit(gaussian_lin_noise_f, Peak511_x_data, Peak511_y_data, Peak511_σ_data.^-2, p0_lin);
Peak511_lin_P = gaussian_lin_511_fit.param .± stderror(gaussian_lin_511_fit)
Peak511_lin_f(x) = gaussian_lin_noise_f(x, gaussian_lin_511_fit.param);
Peak511_lin_quality_R  = quality_of_fit_f(Peak511_lin_f, Peak511_x_data, Peak511_y_data,
                                      Peak511_σ_data, length(gaussian_lin_511_fit.param))

gaussian_const_511_fit = curve_fit(gaussian_const_noise_f, Peak511_x_data, Peak511_y_data, Peak511_σ_data.^-2, p0_const);
Peak511_const_P = gaussian_const_511_fit.param .± stderror(gaussian_const_511_fit)
Peak511_const_f(x) = gaussian_const_noise_f(x, gaussian_const_511_fit.param);
Peak511_const_quality_R  = quality_of_fit_f(Peak511_const_f, Peak511_x_data, Peak511_y_data,
                                      Peak511_σ_data, length(gaussian_const_511_fit.param))


begin
    p_511= plot(spectrum_data[!, :Channel], spectrum_data[!, :Count],
        seriestype=:step,
        xlims=(250, 350), size=(800,500),
        title="511keV peak", label=:none, legend=:topleft
    )
    plot!(250:350, Peak511_quad_f, label="Quadractic: χ²/ndof = $(round(Peak511_quad_quality_R[end], digits=2))")
    plot!(250:350, Peak511_lin_f, label="Linear: χ²/ndof = $(round(Peak511_lin_quality_R[end], digits=2))")
    plot!(250:350, Peak511_const_f, label="Constant: χ²/ndof = $(round(Peak511_const_quality_R[end], digits=2))")
    plot!(250:350, Peak511_f, label="No removal: χ²/ndof = $(round(Peak511_quality_R[end], digits=2))")
end



## 511 edge
p_compton_continuum_511 = plot(spectrum_data[!, :Channel], spectrum_data[!, :Count],
    seriestype=:step,
    xlims=(0, 275),
    xticks=0:50:275,
    title="Compton continuum and edge of 511keV"
)



## Energy calibration




## Testing linear fit
