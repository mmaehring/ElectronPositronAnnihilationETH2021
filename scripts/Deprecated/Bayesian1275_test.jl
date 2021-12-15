
##
using StatsPlots
using BAT
import IntervalSets
using ValueShapes
using CSV
using DataFrames
using StatsBase
using Distributions

path = "data\\3.3\\3.3.2_G\\3.3.2.txt"; # data is stored in 3.3.2 -> this was the acquisition
spectrum_data = CSV.read(path, DataFrame);
Peak1275_data = spectrum_data[725:825, :];

Peak1275_hist = StatsBase.fit(Histogram, Peak1275_data[!,:Channel],
                              weights(Peak1275_data[!,:Count]), nbins=length(Peak1275_data[!,:Count])
)

## Bayesian
function fit_function(p::NamedTuple{(:a, :mu, :sigma)}, x::Real)
    p.a[1] * pdf(Normal(p.mu[1], p.sigma), x)
end

likelihood = let h = Peak1275_hist, f = fit_function
    # Histogram counts for each bin as an array:
    observed_counts = h.weights

    # Histogram binning:
    bin_edges = h.edges[1]
    bin_edges_left = bin_edges[1:end-1]
    bin_edges_right = bin_edges[2:end]
    bin_widths = bin_edges_right - bin_edges_left
    bin_centers = (bin_edges_right + bin_edges_left) / 2

    params -> begin
        # Log-likelihood for a single bin:
        function bin_log_likelihood(i)
            # Simple mid-point rule integration of fit function `f` over bin:
            expected_counts = bin_widths[i] * f(params, bin_centers[i])
            logpdf(Poisson(expected_counts), observed_counts[i])
        end

        # Sum log-likelihood over bins:
        idxs = eachindex(observed_counts)
        ll_value = bin_log_likelihood(idxs[1])
        for i in idxs[2:end]
            ll_value += bin_log_likelihood(i)
        end

        # Wrap `ll_value` in `LogDVal` so BAT knows it's a log density-value.
        return LogDVal(ll_value)
    end
end

prior = NamedTupleDist(
    a = Weibull(30, 1800),
    mu = [IntervalSets.ClosedInterval(750, 800)],
    sigma = Weibull(10, 25)
);

posterior = PosteriorDensity(likelihood, prior);

samples = bat_sample(posterior,
                    MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^4, nchains = 4)
).result;

plot(
    samples, :mu,
    mean = true, std = true, globalmode = true, marginalmode = true,
    nbins = 50, title = "Marginalized Distribution for mu[1]",
    size=(800,500)
)
plot(
    samples, (:(mu), :sigma),
    mean = true, std = true, globalmode = true, marginalmode = true,
    nbins = 50, title = "Marginalized Distribution for μ and σ", size=(800,500)
)

samples_mode = mode(samples)
findmode_result = bat_findmode(posterior, MaxDensityNelderMead(init = ExplicitInit([samples_mode])));
fit_par_values = findmode_result.result[]

plot(
    Peak1275_hist,
    color=1, linewidth=2, fillalpha=0.0,
    st = :steps, fill=false, label = "Data",
    title = "Data, True Model and Best Fit"
)

gaussian_f_new(x) = fit_par_values[1] * exp( -1*(x - fit_par_values[2][1])^2 / ( 2*fit_par_values[3]^2) )

plot!(724.0:1.0:825.0, x -> gaussian_f_new(x), label = "Truth")

# plot!(724.0:0.01:825.0, fit_function, samples)
