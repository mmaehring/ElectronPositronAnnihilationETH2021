#=
In this file i assume we have measured the magnitude of the restoring force
from a simple spring.
=#
## Read in data
using DataFrames
using CSV

data_path_ = "C:\\Users\\marcu\\Pwogwamming\\JuliaThings\\ETH\\some_error_variables.csv"
ydata = CSV.File(data_path,
                        delim=',',
                        ignorerepeated=false,
                        missingstring="NA") |> DataFrame

## Begin processing data
using Measurements
using Unitful

xdata = [i ± 0.25 for i in 1:1:10] # cm
ydata[!, "Combined"] = measurement.(string.(ydata[!, "Value"]) .* "+-" .* string.(ydata[!, "Error"])) # N
ydata = ydata[!, :Combined] |> Vector

## Take a first look at the data, is something obviously wrong?
using StatsPlots
values_and_fits = plot(xdata, ydata,
                        markersize = 4,
                        seriestype=:scatter,
                        # xlabel = "Displacement [cm]",
                        xticks = 1:1:10,
                        ylabel = "Restoring Force [N]",
                        label = "Measured Data",
                        legend = :bottomright
)

## First fit the data with a Least Squares fit (with errors??)
using LsqFit
model = @. model(x, p) = p[1] * x + p[2] # Linear Model -> we expect p[2] = 0

function jacobian_model(x, p) # implement jacobian of linear function
    J = Array{Float64}(undef, length(x), length(p))
    @. J[:,1] = x     #dmodel but derive by p[1]
    @. J[:,2] = 1     #dmodel but derive by dp[2],
    return J
end

# Initial values
p0 = [0.5, 0.5]

# Perform the fit and read out parameters
# fitted = curve_fit(model, jacobian_model,
#                     Measurements.value.(xdata), Measurements.value.(ydata),
#                     p0
# )
# p1, p2 = fitted.param

# Errors of the ydata
weights_errors = Measurements.uncertainty.(ydata)

# Compute with accounting for error in y
fitted = curve_fit(model, jacobian_model,
                    Measurements.value.(xdata), Measurements.value.(ydata),
                    weights_errors,
                    p0
)
p1err, p2err = fitted.param
residual_vals = fitted.resid .± Measurements.uncertainty.(ydata)

# Plot suggested values
new_x = collect(1:0.1:10)
# new_y(x) = @. p1 * x + p2
# plot!(new_x, new_y,
#     label = "Fitted -> No error"
# )
new_y_err(x) = @. p1err * x + p2err
plot!(new_x, new_y_err,
    label = "Fitted -> Errors"
)

## Plots residuals and pulls along with rest of thing
residual_plot = plot(xdata, residual_vals,
                    seriestype=:scatter,
                    label = "Residuals",
                    xlabel = "Displacement [cm]",
                    xticks = 1:1:10,
                    ylabel = "Fit residuals",
                    legend = :bottomleft,
                    legendfontsize=6
)
plot!((x->0), [1, 10], linestyle = :dash, label=:none, linecolor=:red)


l = @layout [a{0.6h};
             b{0.4h}]
plot(values_and_fits, residual_plot, layout = l)


## Plot pulls and distribution of ...
using StatsBase
using Distributions

l_pull = @layout[a{0.75h}; b]

pulls = residual_vals ./ Measurements.uncertainty.(ydata)
pull_plot = plot(xdata, pulls,
                 seriestype = :scatter,
                 label = "Pulls",
                 xlabel = "Displacement [cm]",
                 ylabel = "Pulls",
                 legend = :top,
                 xticks = 1:1:10,
                 yrange = [-1.5, 2]
                 # xrange = [-1, 1]
)

actual_mean_array = zeros(length(pulls)) .+ mean(pulls)
# plot!(1:1:10, Measurements.value.(actual_mean_array), label = "Actual Average")
plot!(1:1:10, Measurements.value.(actual_mean_array), ribbon = zeros(10) .+ Measurements.uncertainty.(mean(pulls)) ,
      fillalpha = 0.25, c = 1, lw = 2,
      label = "Actual Average w. Error",
      linecolor = :red,
      fillcolor = :red,
      legendfontsize = 6
)
plot!(1:1:10, zeros(10), ribbon = ones(10) , fillalpha = 0.25, c = 1, lw = 2,
      label = "1σ around 0"
)


# pull_distribution = histogram(Measurements.value.(pulls)
# pull_histogram = StatsBase.fit(Histogram, Measurements.value.(pulls), nbins=4)
# pull_histogram.weights

hist = histogram(Measurements.value.(pulls), bins = range(-1, 1; length = 6),
    weights = ones(length(pulls)) / length(pulls),
    xlabel = "Pulls",
    ylabel = "PDF"
)

fitted_distribution = Distributions.fit(Distributions.Normal, Measurements.value.(pulls))

plot!(fitted_distribution)

plot(pull_plot, hist, layout = l_pull)


## Hypothesis testing the fit

# H₀ : There is no linear relationship with the specified parameters between the data points.
χ² = sum( (residual_vals ./  Measurements.uncertainty.(residual_vals)).^2 )
degrees_of_freedom = length(ydata) - 2 # - 2 for the fit parameters

χ²_dof = χ² / degrees_of_freedom #


χ²_distr = Distributions.Chisq(degrees_of_freedom - 1)

α = 0.05
χ²_cutoff = Distributions.quantile(χ²_distr, 1 - α) # cutoff for

if χ² < χ²_cutoff
    print("Reject H₀! The fit is good at a 0.05 significance. ")
end

# Best p - value for this χ²:
p = 1 - Distributions.cdf(χ²_distr, Measurements.value(χ²))
