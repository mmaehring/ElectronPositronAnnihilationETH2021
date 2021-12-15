## Marcus Mähring 2021

## Loading packages
begin
    using DrWatson
    # Plotting
    using StatsPlots, UnitfulRecipes
    # Dealing with data
    using Unitful, Measurements, DataFrames, CSV, Statistics, StatsBase
    # Fitting
    using LsqFit, Distributions
    # additional utilities
    using SpecialFunctions, Interpolations
end

## Helper functions
function quality_of_fit_f(fit_function, x, y, σ, dependent_vars = 3)
        χ² = sum((fit_function.(x) .- y).^2 ./ σ.^2)
        dof = length(x) - dependent_vars
        χ²dof = χ²/dof
        return dof, χ², χ²dof
end

function ratio_test(data, model)
    return 0;
end

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
        title="Entire spectrum - Uncalibrated", ylabel="Detected counts",
        size=(1200,800), label="Counts", dpi=900
    )
    fontsize=7
    annotate!(80, 2.07e4, text("Backscattering", :green, fontsize))
    annotate!(215, 4.6e3, text("511 keV Compton edge", :black, fontsize, rotation = -70.0))
    annotate!(120, 1.1e4, text("Photoelectric effect", :red, fontsize))
    annotate!(310, 2.175e4, text("511 keV peak", :darkblue, fontsize))
    annotate!(600, 1.7e3, text("1275 keV Compton edge", :purple, fontsize))
    annotate!(780, 2.5e3, text("1275 keV peak", :orange, fontsize))

    # savefig("Spectrum_Big.pdf")
end

p_uncal_log = plot(spectrum_data[!, :Channel], log.(10, spectrum_data[!, :Count]),
    seriestype=:step,
    xlims=(0,950),
    xticks=0:100:900,
    title="Entire spectrum, log",
    yticks=0:0.5:5, label="Counts"
)

## 1275 peak   # Frequentist  -> there is very little noise here -> no background fit

p0 = [2250.0, 780.0, 20.0]

Peak1275_data = spectrum_data[725:825, :];
Peak1275_x_data = Peak1275_data[!, :Channel]
Peak1275_y_data = Peak1275_data[!, :Count]

Peak1275_fit = curve_fit(gaussian_f, Peak1275_x_data, Peak1275_y_data, p0)
Peak1275_f(x) = gaussian_f(x , Peak1275_fit.param)

values_1275 = Peak1275_fit.param
sigma1275 = stderror(Peak1275_fit); # cov_mat1275 = estimate_covar(Peak1275_fit);
margin_of_error1275 = margin_error(Peak1275_fit, 0.05);
confidence_inter1275 = confidence_interval(Peak1275_fit, 0.05);

p_1275= plot(spectrum_data[!, :Channel][720:825], spectrum_data[!, :Count][720:825],
            seriestype=:step, title="1275 keV peak", label="Recorded data"
)
plot!(720:825, Peak1275_f, label="Fit, μ = $(round(values_1275[2], digits=2)) ± $(round.(sigma1275[2], digits=2))")

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
