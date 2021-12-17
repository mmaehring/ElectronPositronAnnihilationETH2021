
## Loading packages
begin
    using LsqFit
    using Distributions
    using Statistics
    using StatsPlots
end

## Importing data
data = [2144, 2187, 2237, 2256, 2345, 2182, 2225, 2203, 2302, 2251, 2132, 2236,
        2212, 2209, 2283, 2144, 2191, 2211, 2258, 2239, 2195, 2149, 2137, 2157];
adj_data = [2230, 2202, 2211, 2294, 2162, 2232, 2234, 2225, 2275,
            2238, 2154, 2250, 2148, 2189, 2183, 2218, 2221];

## Fitting distribution

pdf_f = Distributions.fit_mle(Poisson, data)
adj_pdf_f = Distributions.fit_mle(Poisson, adj_data)

## Plotting distributions
begin
    p1 = plot(pdf_f, label="First data, λ=$(round(params(pdf_f)[1], digits=2))", xlims=(2000, 2400));
    plot!(adj_pdf_f, label = "Second data λ=$(round(params(adj_pdf_f)[1], digits=2))",
    legend=:topleft, ylabel="Poisson probability");
end

gaussian_ratio = [pdf(pdf_f, i) / pdf(adj_pdf_f, i) for i in 2000:1:2400];
gaussian_ratio_inv = gaussian_ratio.^-1;

begin
    p2 = scatter(2000:5:2400, gaussian_ratio[begin:5:end], label=:none,
    xlabel="Detected Counts", yticks=0.6:0.2:1.4, markersize=3,
    ylabel="Quotient blue/red", color=:purple
    )
    hline!([1], color=:black, label=:none)
end

begin
    l = @layout [a{0.75h};
    b{0.25h}]

    plot(p1, p2, layout=l, dpi=900, size=(700, 450))
end

# savefig("plots\\PoissonQuotient900dpi.pdf")
