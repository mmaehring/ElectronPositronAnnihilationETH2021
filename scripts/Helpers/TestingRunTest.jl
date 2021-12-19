
## Packages
begin
    using Statistics
    using StatsPlots
    using Measurements
    using LsqFit
    using Distributions
    using Logging
    using Random
end

## Models
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
function ratio_test(data, model_data)
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

## Data
begin
    x = collect(1.0:1.0:12.0)
    y_real = 0.3*x.^2 .+ 1.25
    yerr = 5.5
    yerr = ones(length(y_real)) * yerr
end

begin
    p0 = [12.0, 5.0];
    model(x,p) = @. p[1]*x + p[2]
    fit_model = curve_fit(model, x, y_real, 1 ./yerr.^2, p0)

    yfunc(x) = fit_model.param[1]*x .+ fit_model.param[2];
    chisq = χ²_calc(yfunc, x, y_real, yerr)[2];
    dofχ = length(y) - 3;
    println(χ²_test(chisq, dofχ))
end

begin
    plot(x, yfunc(x), ls=:dashdot, legend=:topleft, label="Fitted line")
    scatter!(x, y_real .± yerr, label="Real data")
    annotate!(5, 35, text("Number of runs: $(compute_runs(y_real, yfunc(x)))"))
end


## Computing runs
begin # polynomial
    Random.seed!(314)
    x = collect(1.0:0.5:12.0)
    y_real = 0.3*x.^2 .+ 1.25 .+ 2*randn(length(x))
    yerr = 5.5
    yerr = ones(length(y_real)) .* yerr

    p0 = [12.0, 5.0];
    model(x,p) = @. p[1]*x + p[2]
    fit_model = curve_fit(model, x, y_real, 1 ./yerr.^2, p0)
    yfunc(x) = fit_model.param[1]*x .+ fit_model.param[2];
    plot(x, yfunc(x), ls=:dashdot, legend=:topleft, label="Fitted line")
    scatter!(x, y_real .± yerr, label="Real data")
    annotate!(5, 35, text("Number of runs: $(compute_runs(y_real, yfunc(x)))", 10))
    annotate!(5, 25, text("Probability: $(round(prob_r_runs(y_real, yfunc(x)), digits=3)*1e2)%", 10))
    r, rexp, V, σ, N, Na = ratio_test(y_real, yfunc(x))
    # @info "$(100 - 100*cdf(Distributions.Normal(), (rexp - r)/σ)) %"
    annotate!(7.5, -3, text("We observe $(r) runs but expect $(rexp±σ) => $(round((rexp-r)/σ, digits=2)) σ", 9))
    annotate!(7.5, -5, text("Good at a one tailed test with probability $( round(100 - 100*cdf(Distributions.Normal(), (rexp - r)/σ), digits=2)) %", 9))
end
