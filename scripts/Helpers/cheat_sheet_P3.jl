
## DEALING WITH PHYSICAL QUANTITIES & ERROR PROPAGATION
using Measurements
# https://juliaphysics.github.io/Measurements.jl/stable/
# The method used to handle functional correlation is described in this paper:
# M. Giordano, 2016, "Uncertainty propagation with functionally correlated quantities",
# arXiv:1610.08716 (Bibcode: 2016arXiv161008716G)


# explicitly declare a measurement
a = measurement(4.5, 0.1)

# "Measurement("...")" creates object from a string.

# declare it in without explicit declaration, using typecasting with \pm + tab
b = 3.8 ± 0.4

# the errors are propagated via ... [reference to documentation -> method]
c = 2a + b

# Note how it also recognizes how a measurement should work with itself (Functional correlation):
x = 8.4 ± 0.7
x - x == 0
x / x == 1
x*x*x - x^3 == 0

# BUT not always exactly
sin(x)/cos(x) - tan(x)

# Works with non-trivial expressions
x = 8.4 ± 0.7;
y = 9.36 ± 1.02;
log(2x^2 - 3.4y)

# Access elements of a value with errors
value = Measurements.value(x)
uncertainty = Measurements.uncertainty(x)
measurement("123.4e-1 +- 0.056e1")
# Also:
# Support for most mathematical operations available in Julia standard library and
# special functions from SpecialFunctions.jl package, involving real and complex numbers.

#= the @uncertain macro
can be used to propagate uncertainty in arbitrary real or complex functions of
real arguments, including functions not natively supported by this package. =#
using SpecialFunctions
@uncertain (x -> complex(zeta(x), exp(eta(x)^2)))(2 ± 0.13)

@uncertain log(9.4 ± 1.3, 58.8 ± 3.7)

## DataFrames
using DataFrames
# https://dataframes.juliadata.org/stable/man/basics/

DataFrame() # empty

DataFrame(A=1:3, B=5:7, fixed=1)

DataFrame("customer age" => [15, 20, 25],
                 "first name" => ["Rohit", "Rahul", "Akshat"])

mat = [1 2 4 5; 15 58 69 41; 23 21 26 69];
nms = ["a", "b", "c", "d"];
DataFrame(mat, nms)

### READING FROM CSV FILES
using CSV
# PATH = "..."
# german_ref = CSV.read(PATH, DataFrame)

## SYSTEMS OF LINEAR EQUATIONS (W UNCERTAINTIES)
using LinearAlgebra # (+ Measurements.jl)
A = [(14 ± 0.1) (23 ± 0.2); (-12 ± 0.3) (24 ± 0.4)] # Matrix with uncertainties (why not?)
b = [(7 ± 0.5), (-13 ± 0.6)] # Vector with uncertainties

# Find solution to A*x = b -> propagates error correctly
x = A\b
# Vector algebra also works with uncertainties
det(A)
dot(x, b)

## UNITS
using Unitful
# Has the standard SI units

# Distinguish units
1u"kg" == 1000u"g"
1u"kg" === 1000u"g"
1u"kg" === 1u"kg"

# You can define your own units
const °C = u"°C"
const °F = u"°F"

# Has unit conversion
uconvert(°C, 212°F)
1u"m" |> u"cm" # convert meters to centimeters
1u"eV" |> u"J" # electron volts to joules
# 1u"eV" |> u"kg" # DimensionError: kg and eV are not dimensionally compatible

# using Unitful.DefaultSymbols imports
# 𝐋,𝐌,𝐓,𝐈,𝚯,𝐉,𝐍 for length, mass, time, current, temperature, luminosity, and amount, respectively.

#=
Unitful.promote_to_derived()
Defines promotion rules to use derived SI units in promotion for common dimensions of quantities:
  J (joule) for energy
  N (newton) for force
  W (watt) for power
  Pa (pascal) for pressure
  C (coulomb) for charge
  V (volt) for voltage
  Ω (ohm) for resistance
  F (farad) for capacitance
  H (henry) for inductance
  Wb (weber) for magnetic flux
  T (tesla) for B-field
  J*s (joule-second) for action
=#

### Interplay with Measurements.jl
# Ohm's Law
(50 ± 1)*u"Ω" * (13 ± 2.4)*1e-2*u"A" |> u"V"
# (I've added the explicit conversion to Volts as we otherwise get Ω⋅R)
2pi*sqrt((5.4 ± 0.3)*u"m" / ((9.81 ± 0.01)*u"m/s^2")) # Pendulum's period

# Remove Units
1u"kN*m"/4u"J" |> NoUnits

## CONSTANTS
using PhysicalConstants
import PhysicalConstants.CODATA2018: c_0, ε_0, μ_0, h, ħ
# Vacuum permativity
ε_0 # more info
2*ε_0 # just the value

# Vacuum permeability
μ_0

# Get more information about the constants (unvertainties, relative)
PhysicalConstants.CODATA2018.NewtonianConstantOfGravitation
PhysicalConstants.CODATA2018.SpeedOfLightInVacuum

# Example, ħ and Measurements.jl
measurement(BigFloat, ħ)
measurement(BigFloat, ħ) / (measurement(BigFloat, h) / (2 * big(pi)))

## RANDOM NUMBERS
using Random
# Two random integers (they get big)
rand(Int, 2)

# Random value between 2 and 3
rand((2, 3))

# Two random values between 2 and 3
rand((2, 3), 2)

# 2x2 Matrix of values between 2 and 3
rand((2, 3), (2,2))

rand(Float64, (2, 3)) # 2x3 matrix with values between 0 and 1

# MersenneTwister: The Mersenne Twister is a pseudorandom number generator (PRNG)
# bitrand(MersenneTwister(1234), 10) # 10 element bit vector with seed 1234

## STATISTICS FUNCTIONS
using Statistics
# STATS and uncertainties in arrays (req.: Measurements.jl)
A = [1.13 ± 0.14, 2.98 ± 0.45, 5.46 ± 0.97]
sum(A) # Sum of values
mean(A) # The arithmetic mean
std(A) # The standard deviation
var(A) # The variance
cor(A) # The Pearson correlation matrix of a matrix X
cov(A) # The Covariance of a matrix (Cov(x), x vector, gives variance and divides it with n-1, n=length(x))
median(A)

## FITTING
using LsqFit # Levenberg-Marquardt algorithm for non-linear fitting.
# https://julianlsolvers.github.io/LsqFit.jl/latest/

# DEINE a two-parameter exponential model
# x: array of independent variables
# p: array of model parameters
# model(x, p) will accept the full data set as the first argument `x`.
# This means that we need to write our model function so it applies
# the model to the full dataset. We use `@.` to apply the calculations
# across all rows.
@. model(x, p) = p[1]*exp(-x*p[2])

# some example data
# xdata: independent variables
# ydata: dependent variable
xdata = range(0, stop=10, length=20)
ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))

#Before fitting the data, we also need a initial value of parameters for curve_fit().
p0 = [0.5, 0.5]

fit = curve_fit(model, xdata, ydata, p0)
param = fit.param # equivalent to coef(fit)
# fit is a composite type (LsqFitResult), with some interesting values:
#	dof(fit): degrees of freedom
#	coef(fit): best fit parameters
#	fit.resid: residuals = vector of residuals
#	fit.jacobian: estimated Jacobian at solution

# If we know ranges for the variables (also: greater than zero? [0, Inf])
# Optional upper and/or lower bounds on the free parameters can be passed as an argument.
# Bounded and unbouded variables can be mixed by setting `-Inf` if no lower bounds
# is to be enforced for that variable and similarly for `+Inf`
lb = [1.1, -0.5]
ub = [1.9, Inf]
p0_bounds = [1.2, 1.2] # we have to assign initial values inside the bounds

fit_bounds = curve_fit(model, xdata, ydata, p0_bounds, lower=lb, upper=ub)

###  Errors on fits - estimations

# We can estimate errors on the fit parameters,
# to get standard error of each parameter:
sigma = stderror(fit)
# also the covariance matrix
cov_mat = estimate_covar(fit)
# to get margin of error and confidence interval of each parameter at 5% significance level:
margin_of_error = margin_error(fit, 0.05)
confidence_inter = confidence_interval(fit, 0.05)

#= To get the confidence interval at 10% significance level, run
confidence_interval(fit, alpha)
which essentially computes:
estimate parameter value ± (standard error * critical value from t-distribution). =#
confidence_interval_value_10 = confidence_interval(fit, 0.1)  # (lower and upper bound for each param)

# The finite difference method is used above to approximate the Jacobian.
# Alternatively, a function which calculates it exactly can be supplied instead.
function jacobian_model(x,p)
    J = Array{Float64}(undef, length(x), length(p))
    @. J[:,1] = exp(-x*p[2])     #dmodel/dp[1]
    @. @views J[:,2] = -x*p[1]*J[:,1] #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
    return J
end
fit = curve_fit(model, jacobian_model, xdata, ydata, p0)

# Multivariate regression
#=There's nothing inherently different if there are more than one variable entering the problem.
We just need to specify the columns appropriately in our model specification: =#
@. multimodel(x, p) = p[1]*exp(-x[:, 1]*p[2]+x[:, 2]*p[3])


### To perform WEIGHTED LSTQ we can add the weight parameter (wt) to perform weighted LSTQ
# This seems kind of complicated. Should read up on more.
# https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/#Weighted-Least-Squares-1
# https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/#Estimate-the-Optimal-Weight-1

#= NOTE
curve_fit() assumes the weight as the inverse of the error covariance matrix rather than the ratio
of error covariance matrix, i.e. the covariance of the estimated parameter is calculated as covar = inv(J'*fit.wt*J). =#

# fit = curve_fit(m, tdata, ydata, wt, p0) -> wt be vector of y errs (basic)
# fit = curve_fit(m, tdata, ydata, wt, p0) -> wt be covariance matrix of values



## MINIMIZATION
using Optim
# https://github.com/JuliaNLSolvers/Optim.jl

# Define the rosenbrock function in two variables. (classic function for testing stiff solvers)
rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
# Minimize the function with initial parameter values (0,0) and using BFGS
# which does not need matrix inversion and thus has an algorithmic complexity of 𝒪(n²) instead of 𝒪(n³) for Newtons method.
# If one has the analytic hessian and the dimension of the problem is low, then second order methods such as Newtons is the method of choice.
# https://en.wikipedia.org/wiki/Broyden-Fletcher-Goldfarb-Shanno_algorithm
result1 = optimize(rosenbrock, zeros(2), BFGS())
result2 = optimize(rosenbrock, zeros(2), Newton())


# Minimization with errors / units?? ↯ DOES NOT WORK WITH THE METHODS BELOW
# rb_1 = 1.0 ± 0.05
# rb_2 = 100.0 ± 13.5
# rosenbrock_errors(x) = (rb_1 - x[1])^2 + rb_2 * (x[2] - x[1]^2)^2
# result3 = optimize(rosenbrock_errors, zeros(2), BFGS())

# rb_1 = 1.0u"kg"
# rb_2 = 100.0u"kg"
# rosenbrock_units(x) = (rb_1 - x[1])^2 + rb_2 * (x[2] - x[1]^2)^2
# result4 = optimize(rosenbrock_errors, [0u"kg", 0u"kg"], BFGS())

## PLOTTING
# Basic plots
using Plots
x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
y = 3x .+ 0.5
z = x -> x^2
plot(x, y)
plot(z)
plot(x -> x^3, -4, 4)
plot(rand(10))     # 1 series... x = 1:10
plot([sin,cos], 0:0.1:π)
plot(sin,
    (x->begin sin(2x) end), # equivalent to x -> sin(2x) for some reason
    0, 2π,
    line = 4, leg = false,
    fill = (0, :orange)
)

# plot creates new canvas
# plot! modifies canvas

#

# Bar plot with options
#=
Plots.bar(
  df[!, dice][1:faces],
  label="Recorded Data",
  xlabel="Result",
  ylabel="Occurances",
  show = true,
  title="Occurance of parameters for the " * dice
) =#

# Plotting some uncertainties (both in x and y directions)  # OBS : requires Measurements.jl
plot(sin, [x ± 0.1 for x in 1:0.2:10], size = (1200, 800))

# Plotting some "measured" values, fitting them to a plot
# We have: y = p1 * x + p2 as the theorized model
xdata = [x ± 0.05 for x in 0:0.5:5]
ydata = ([2.75*i + 0.5 for i in xdata] + 0.25*randn(length(xdata)))

plot(xdata, ydata,
    xlabel="xdata",
    ylabel="ydata",
    seriestype=:scatter,
    markersize=3.25,
    dpi=300,
    size=((600, 400)),
    title=(nothing)
)


model = @. model(x, p) = p[1]*x + p[2]

# the jacobian for the model above -> reduces computation time and accuary as we don't need to use forward differentiation.
function jacobian_model(x,p)
    J = Array{Float64}(undef, length(x), length(p))
    @. J[:,1] = x     #dmodel/dp[1]
    @. @views J[:,2] = 1 #dmodel/dp[2], thanks to @views we don't allocate memory for the J[:,1] slice
    return J
end

lb  = [0, -Inf]  # define the lower boundary

fit = curve_fit(model, jacobian_model,
                Measurements.value.(xdata),
                Measurements.value.(ydata),
                p0,
                lower = lb
) # perform fit -> need to figure out how to use errors (weighted)

p1, p2 = fit.param                                                              # get the fitted values
new_and_improved_model(x) = p1 * x + p2                                         # construct function with the fitted values
new_linspace = collect(0:0.01:5)                                                # create new linspace
plot!(new_linspace, new_and_improved_model)                                     # plot the approximated function in the same plot as the data

residues = residuals(fit)
plot(xdata, residues, seriestype=:scatter, title="Residuals")
plot!([0:0.5:5], ones(length(xdata)) .- 1, )


### Advanced Plots : Plotting residue and pull of a fit under
l = @layout [a{0.6h};
             b{0.4h}]
plot1 = plot(xdata, ydata,
    ylabel="ydata [unit]",
    seriestype=:scatter,
    markersize=3.25,
    dpi=250,
    size=((1200, 800)),
    title = "Data",
    label = "Measured data",
    legend = :bottomright,
    xticks = 0:0.5:5,
    legendfontsize=7
)
plot!(new_linspace, new_and_improved_model.(new_linspace), label="Fit")


### MORE COLOR SCHEMES
# import ColorSchemes
# https://juliagraphics.github.io/ColorSchemes.jl/stable/
plot2 = plot(xdata, residues,
            seriestype=:scatter,
            ylabel="Residuals",
            xlabel = "xdata [unit]",
            label = "Residuals",
            # palette = ColorSchemes.spring.colors
            xticks = 0:0.5:5,
            legendfontsize=6,
            legend = :topright,
)
plot!([0:0.5:5], ones(length(xdata)) .- 1, label = "0")

plot(plot1, plot2, layout = l)

# More advanced statistics plotting tools - > provides plotting recipes
# using StatsPlots # (also imports Plots.jl)
#= Don't know if this is actually from StatsPlots
Most importantly this package provides:
 ⊗ histogram/histrogram2d
 ⊗ Kernel density estimation plots (marginalkde)
=#

### LOGPLOTS - NOTE IMPORTANT -> don't include 0 in the linspace you want to plot!
f(x) = x^-1.7
plot(f, 1e-3, 3, scale=:log10, title="log-log")
plot(f, 1e-3, 3, xscale=:log10, title="log-log")
plot(f, 1e-3, 3, yscale=:log10, title="log-log")

### RIBBONS (FILLING BETWEEN POINTS)
x = collect(range(0, 2, length= 100))
y1 = exp.(x)
y2 = exp.(1.3 .* x)

plot(x, y1, fillrange = y2, fillalpha = 0.35, c = 1, label = "Confidence band", legend = :topleft)

mid = (y1 .+ y2) ./ 2   #the midpoints (usually representing mean values)
w = (y2 .- y1) ./ 2     #the vertical deviation around the means

plot(x, mid, ribbon = w , fillalpha = 0.35, c = 1, lw = 2, legend = :topleft, label = "Mean")
plot!(x,y1, line = :scatter, msw = 0, ms = 2.5, label = "Lower bound")
plot!(x,y2, line = :scatter, msw = 0, ms = 2.5, label = "Upper bound")


## DISTRIBUTIONS
using Distributions
using StatsPlots
plot(Normal(3,5), fill=(0, .5,:orange)) # -> plot recipe from StatsPlots , doesn't work for Plots.jl


d = Normal() # Create a normal distribution
x = rand(d, 100) # Sample the distribution 100 times

### Get pdf, cdf, quantile

# Median (50th percentile) and the 95th percentile for the standard-normal distribution are given by:
quantile.(Normal(), [0.5, 0.95])

# To draw random samples from a normal distribution with mean 1 and standard deviation 2, write:
rand(Normal(1, 2), 100)

### Some distributions
# Binomial(p) # Discrete univariate
# Cauchy(u, b)  # Continuous univariate
# Multinomial(n, p) # Discrete multivariate
# Wishart(nu, S) # Continuous matrix-variate

# get fieldnames:
fieldnames(Cauchy)

# ESTIMATE PARAMETERS
# Here we fit randomly sampled data (100 points) from the Normal distribution d above to a Gaussian
x = rand(d, 100)
Distributions.fit(Normal, x)
# true values [0.0, 1.0]


## HYPOTHESIS TESTING
degrees_of_freedom = 4 # (number of possible sixes rolling 3 dice - 0, 1, 2, or 3 sixes)
# probability distribution describing [?????]
χ²_distr = Distributions.Chisq(degrees_of_freedom - 1)
# we subtract 1 to get a more accurate (and more restrictive) estimate for the dofs

# the χ² of the data points
# χ²_val = Σ (observed - expected)² / expected
#
# general form- for regressions, see below.

χ²_val = 23.367
α = 0.05

#= If we are interested in a significance of α (here 0.05 = 5%) for a χ² distribution with a certain
amount of degrees of freedom; then

IF  χ² ≥ Distributions.quantile(χ²_distr, 1-α) = Z(α) = 7.815
 with Z(α) = (inverse survival function).

THEN we can reject the H₀ hypothesis.

EG: in this example, the H₀ hypothesis is: "The dice are fair"
which we subsequently reject considering χ²_val = 23.367  ≥  7.815
The H₀ hypothesis is the dice are fair, so to the answer to "are the dice fair?" is NO
=#

cutoff = Distributions.quantile(χ²_distr, 1 - α)

if χ²_val ≥ Distributions.quantile(χ²_distr, 1 - α)
    print("Reject H₀ ⟹ The dice are not fair \n")
end

# P-value of the χ² value:
p_val_χ² = (1 - Distributions.cdf(χ²_distr, χ²_val))
# very small : this makes sense, as H₁ doesn't really describe what we see at all.

### NOTE : Let us do the same thing for a fit of data to a suspected (linear) function.
xdata = [x ± 0.05 for x in 0:0.5:5]
ydata = [2.75*i + 0.5 for i in xdata] + 0.25 * randn(length(xdata))
model = @. model(x, p) = p[1]*x + p[2]
fit = curve_fit(model,
                Measurements.value.(xdata),
                Measurements.value.(ydata),
                p0
) # simplified fitting -> ignore bounds and numerical jacobian. (complete: see above)
# also: identical errors so no weighting necessary
p1, p2 = coef(fit)

# TODO
### AND -> residues -> χ² = ∑ᵢ (yᵒᵢ - yᵗᵢ)² / σᵢ²
# here i runs over then bins, yᵒ  is the observed values and yᵗ the theoretical.


## BINNING & SOME MORE HYPOTHESIS TESTING

### BINNING
using StatsBase
pull_plot = StatsBase.fit(Histogram, Measurements.value.(pulls), nbins=6)
histo.weights

### χ² statistics for 200 students SAT scores -> H₀ : They follow a normal distribution
using CSV
path_to_data = "C:\\Users\\marcu\\Pwogwamming\\JuliaThings\\ETH\\test_scores.csv"
test_scores = CSV.File(path_to_data,
                        delim=' ',
                        ignorerepeated=false,
                        missingstring="NA") |> DataFrame
test_scores_array = (test_scores |> Array)
test_scores_vector = reshape(test_scores_array, (200))

# binning : (< -2.0) , (-2.0, -1.5) , (-1.5, -1.0) , (-1.0, -0.5), (-0.5, 0.0), (0.0, 0.5), (0.5, 1.0) , (1.0, 1.5) , (1.5, 2.0), (> 2.0)
# exptected scores -> from probabilities => stems from binning i guess
probabilities = [0.023, 0.044, 0.092, 0.150, 0.191, 0.191, 0.150, 0.092, 0.044, 0.023]
expected_people = [4.6, 8.8, 18.4, 30.0, 38.2, 38.2, 30.0, 18.4, 8.8, 4.6]
measured = [6, 6, 18, 33, 38, 38, 28, 21, 9, 3]

# chisq = sum((measured - expected)²  / expected)
χ² = sum( (measured - expected_people).^2 ./ expected_people)

#= Since the data are divided into 10 bins and we have estimated two parameters,
the calculated value may be tested against the chi-square distribution with 10 -1 -2 = 7 degrees of freedom. =#
α = 0.05
degrees_of_freedom = 10 - 2
χ²_distr = Distributions.Chisq(degrees_of_freedom - 1)
χ²_cutoff = Distributions.quantile(χ²_distr, 1 - α)
#=
For this distribution, the critical value for the 0.05 significance level is 14.07.
Since 2.69 < 14.07, we do not reject the null hypothesis that the data are normally distributed.
=#
χ² < χ²_cutoff
# ⟹ accept hypothesis, reject that this function does NOT model the data to an accuracy of α = 0.05

## MAXIMUM LIKELIHOOD ESTIMATION -> "In statistics, maximum likelihood estimation (MLE)
#is a method of estimating the parameters of an assumed probability distribution,
# given some observed data."

### Example 1 : Light Bulbs
#= Assume we can model the lifespan of a brand of lightbulbs by an exponential distribution
with an unknown parameter λ. Assume we measured lifetimes of 2, 3, 1, 3, 4 years respectively.
Now, let Xᵢ be the lifetime of the ith bulb, and let xᵢ be the value Xᵢ takes. Then each
Xᵢ has a PDF of the form:
f_Xᵢ(xᵢ) = λ exp(-λxᵢ).
Assume the lifetimes of the bulbs are independent, thus the product of the PDFs
gives the joint PDF.
f(x₁, ..., x₅ | λ) = Πᵢ λ ⋅ exp(-λ * xᵢ)
By explicity plugging in values and then differentiating, one finds the optimal λ ≈ 5/13
=#
values = [2, 3, 1, 3, 4]

# LIKELIHOOD,  ℒ
DISTR(λ) = prod(λ .* exp.(-λ .* values)) # pdfs
ℒ(λ) = DISTR(λ)
plot(0:0.01:2, DISTR, xlims=[0, 1.25], ylims=[0, 6e-5])
maximized_likelihood = maximize(DISTR, ones(1), BFGS())
maxim = Optim.maximizer(maximized_likelihood)[1] #  get's the found maximizing value

# Equivalent -> ℓ = log(ℒ(λ)) = Σᵢ log(yᵢ)
ℓ(λ) =  -1 * sum( log.(λ .* exp.(-λ .* values)) )
res = optimize(ℓ, [0.2], BFGS()) # minimize the negative function gives maximum
maxim1 = first(Optim.minimizer(res))

### Normal Distributions
# Assume we have x₁ ... xₙ drawn from a normal distribution N(μ, σ²) with μ, σ unknown.
# Find maximum likelihood estimator for (μ, σ).
#=
Let X1 ... Xn be independent N(μ, σ²) random variables with xᵢ the value the random variable takes.
Then the density for each random variable will be:
f_Xᵢ (xᵢ) = 1/(srqt(2π)σ) ⋅ exp(- (xᵢ - μ)² / (2σ²))
=#
function NormalImpl(μ, σ, x)
    return @. 1/(sqrt(2π) * σ) * exp(-1 * (x - μ)^2 / (2*σ^2))
end

function MLE_Normal(xdata)
    ℒ(x) = prod(NormalImpl(x[1], x[2], xdata))
    maximized_likelihood = maximize(ℒ, ones(2), Newton())
    maxim = Optim.maximizer(maximized_likelihood)
    return maxim
end

function MLE_NLog(xdata)
    ℒ(x) = (sum(log.(NormalImpl(x[1], x[2], xdata))))
    maximized_likelihood = maximize(ℒ, ones(2), Newton())
    maxim = Optim.maximizer(maximized_likelihood)
    return maxim
end

# Here the MLE estimation based on multiplication fails, probably due to small numbers. -> log is the way to go.
some_data = rand(Normal(0, 1), 50) #μ = 0, σ = 1
# MLE_Normal(some_data)
MLE_NLog(some_data)

## Trying different methods to find initial guess for minimization problems:

# Brute-forcing it with a linspace
points = collect(0.01:0.05:5); eval_points = [ℓ(x) for x in points]; # generate sample spacing + evaluate on spacing
min_index = argmin(values)   ; # find minimum index
start_val = points[min_index] # find corresponding x-value
res = optimize(ℓ, [start_val], Newton()) # minimize the negative function gives maximu
maxim1 = first(Optim.minimizer(res))

### Two iterations of Newton and then minimizing -> turning a bad guess into a good one
using FiniteDifferences # package to approximate derivatives
# https://github.com/JuliaDiff/FiniteDifferences.jl
# Exkurs in FiniteDifferences : # which method -> central_fdm, which order -> 5 Which order derivateive? -> 1, which function? -> sin, which point? -> 1
example_first_derivative = central_fdm(5, 1)(sin, 1) - cos(1)
example_second_derivative = central_fdm(5, 2)(sin, 1) + sin(1)

# Newton iterations
# https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
function Newton_iteration(fct, start, iterations)
    iterations += 1
    values = zeros(iterations)
    values[1] = start
    for i in 2:iterations
        approximated_derivative = central_fdm(5, 1)(fct, values[i-1])
        values[i] = abs(values[i-1] - fct(values[i-1]) / approximated_derivative)
    end
    return last(values)
end

### Sample starting point
xₛ = 0.5
iterations = 2
better_guess = Newton_iteration(ℓ, xₛ, iterations)
# Minimize
minimized_value = first(Optim.minimizer(optimize(ℓ, [better_guess], Newton()))) # smaller steps with newton=

## QUADRATURE
using QuadGK # Adaptive gaussian quadrature
# https://github.com/JuliaMath/QuadGK.jl

# General use case:
integral, err = quadgk(x -> exp(-x^2), 0, 1, rtol=1e-8)

# For computing many integrands of similar functions with singularities
x, w = gauss(x -> exp(-x) / sqrt(x), 10, 0, -log(1e-10), rtol=1e-9)

# With Measurements
a = 4.71 ± 0.01;
quadgk(x -> exp(x / a), 1, 7)[1]

# With Unitful
a = (4.71 ± 0.01) * u"V^3"
# Add units to the bounds and your prayers shall be recieved.
quadgk(x -> x^3 / a, 1u"V", 7u"V")[1]

# (don't know how to specify that x is V^4 ⟹ use *1...unit to force it)
# quadgk(x -> x^3 / a, 1, 7)[1] * 1u"V^4"

## NUMERICAL QUADRATURE FROM DISCRETE POINTS
using Trapz

vx=range(0,1,length=100)
vy=range(0,2,length=200)
vz=range(0,3,length=300)

# 3D
M = [x^2+y^2+z^2 for x=vx, y=vy, z=vz] # Define values for function on a 100*200*300 cuboid

I=trapz((vx,vy,vz), M) # Numerically integrate over these discrete points

print("result: ", I)

#=
@trapz range variable expression
Calculates integral of [expression] over [variable] in [range]
Example:
=#
@trapz 0:0.01:1 x x*x

# 1D with errors #TODO

## DIFFERENTIAL EQUATIONS
using DifferentialEquations
# https://diffeq.sciml.ai/stable/#Supporting-and-Citing
f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

#Plot the solution
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) [μm]", label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t), lw=3, ls=:dash, label="True Solution!")

### With units and errors (Unitful / Measurements)
### Boundary value problem
using BoundaryValueDiffEq
const g = 9.81 ± 0.01
L = 1.0
tspan = (0.0,pi/2)
function simplependulum!(du,u,p,t)
    θ  = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L)*sin(θ)
end

# TODO

### A small list of solvers
#=
    AutoTsit5(Rosenbrock23()) handles both stiff and non-stiff equations. This is a good algorithm to use if you know nothing about the equation.
    AutoVern7(Rodas5()) handles both stiff and non-stiff equations in a way that's efficient for high accuracy.
    Tsit5() for standard non-stiff. This is the first algorithm to try in most cases.
    BS3() for fast low accuracy non-stiff.
    Vern7() for high accuracy non-stiff.
    Rodas4() or Rodas5() for small stiff equations with Julia-defined types, events, etc.
    KenCarp4() or TRBDF2() for medium sized (100-2000 ODEs) stiff equations
    RadauIIA5() for really high accuracy stiff equations
    QNDF() for large stiff equations
Source: https://diffeq.sciml.ai/stable/solvers/ode_solve/#ode_solve
=#

### Hamiltonian Problems
using DiffEqPhysics
# https://diffeq.sciml.ai/stable/types/dynamical_types/#dynamical_prob

## MODELINGTOOLKIT
## SOLID STATE PHYSICS
## STATISTICAL PHYSICS
## ANIMATIONS

## ONLINE STATS
using BinningAnalysis

# Binning tools:
# Logarithmic Binning
# Size complexity: O(log(N))
# Time complexity: O(N)
# Full Binning (all bin sizes that work out evenly)

B = LogBinner()
# push!(B, 2.0)
append!(B, test_scores_vector)

B = FullBinner() # <: AbstractVector (lightweight wrapper)
# push!(B, 2.0)
append!(B, [1,2,3])


## Automatic Differentiation
using ReverseDiff
using BenchmarkTools
# https://github.com/JuliaDiff/ReverseDiff.jl
