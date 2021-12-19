begin
    using StatsPlots
    using LsqFit
    using Measurements
end;

## data generation
begin
    x = collect(0:5:25)
    y = 8.321*x .+ 1.5 .+ 0.1*randn(length(x))
    y_err = 10.5
end;
y_err_array = y .± y_err
yarr = ones(length(y))*y_err .* randn(length(y)) # testing different errors

begin
    l = (y.-yarr)
    u = (y.+yarr)
    c_confidence = 1
    plot(x, y, fillrange = u, fillalpha = 0.35, c = c_confidence, label = "Confidence band", linewidth=0)
    plot!(x, y, fillrange = l, fillalpha = 0.35, c = c_confidence, label = :none, linewidth=0)
    plot!(x, y.±y_err, color=:red, legend=:topleft, label="Fitted line")
end


begin
    l = (y.-yarr)
    u = (y.+yarr)
    c_confidence = 1
    plot(x, y.±y_err, color=:red, legend=:topleft, c=:red, label="Fitted line")
    plot(x, (l .+ u) ./ 2, ribbon = (l .- u) ./ 2, fillalpha = 0.35, c = 1, label = "Confidence band")
end
