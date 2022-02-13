

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
    using SpecialFunctions, Interpolations
end

cd(raw"C:\Users\marcu\OneDrive\Desktop\PraktikumIII\e+e-_Annihilation")

## Loading data
begin
    path = "data\\3.4\\3.4.2\\3.4.2.txt"; # data is stored in 3.3.2 -> this was the acquisition
    spectrum_data = CSV.read(path, DataFrame);
    spectrum_data[!, :Error] = sqrt.(spectrum_data[!, :Count])
end

## Plotting spectrum
begin
    x = spectrum_data[!, :ADC]
    y = spectrum_data[!, :Count]
    y_error_array = spectrum_data[!, :Error]
    l = (y .- 2*y_error_array)
    u = (y .+ 2*y_error_array)
    clrs = [1, :red]
    c_confidence = clrs[2]
    α = 0.5
    plot(x, y, fillrange = u, fillalpha = α, c = c_confidence, label = "Confidence band (2σ)", seriestype=:step,
        xlims=(0,400), xticks=0:100:400, xlabel="Channel", linewidth=0,
        title="MCA coincidence spectrum - Uncalibrated", ylabel="Detected counts",
        size=(700,500), dpi=700, titlefontsize=11, xguidefontsize=10, yguidefontsize=10,
        margin = 5mm
    )
    plot!(x, y, fillrange = l, st=:step, fillalpha = α, c = c_confidence, label = :none, linewidth=0)
    plot!(x, y, color=clrs[1], st=:step, legend=:topleft, label="Fitted line", linewidth=0.7)
    # savefig("plots\\3.4.2_MCA_coinc_nocalib.pdf")
end

# Log
begin
    x = spectrum_data[!, :ADC]
    y = spectrum_data[!, :Count]
    y_error_array = spectrum_data[!, :Error]

    y_tmp = y .± y_error_array
    y_log_tmp = log.(10, y_tmp)

    l = Measurements.value.(y_log_tmp) .- 2*Measurements.uncertainty.(y_log_tmp)
    u = Measurements.value.(y_log_tmp) .+ 2*Measurements.uncertainty.(y_log_tmp)
    clrs = [1, :red]
    c_confidence = clrs[2]
    α = 0.5
    plot(x, Measurements.value.(y_log_tmp), fillrange = u, fillalpha = α,
        c = c_confidence, label = "Confidence band (2σ)", seriestype=:step,
        xlims=(0,400), xticks=0:100:400, xlabel="Channel", linewidth=0,
        title="MCA coincidence spectrum, log - Uncalibrated", ylabel="Detected counts, log₁₀",
        size=(700,500), dpi=700, titlefontsize=11, xguidefontsize=10, yguidefontsize=10,
        margin = 5mm, ylims=(1,3.5)
    )
    plot!(x, Measurements.value.(y_log_tmp), fillrange = l, st=:step, fillalpha = α, c = c_confidence, label = :none)
    plot!(x, Measurements.value.(y_log_tmp), color=clrs[1], st=:step, legend=:topleft, label="Fitted line", linewidth=0.7)
    # savefig("plots\\3.4.2_MCA_coinc_nocalib_LOG.pdf")
end



## Converting ADC to energies
E_ADC(x) = (-11.4 ± 3.5) + (0.617 ± 0.0037) * x
ADC_E(x) = (x - (-11.4 ± 3.5)) / (0.617 ± 0.0037)

begin
    x = spectrum_data[!, :ADC]
    y = spectrum_data[!, :Count]
    y_error_array = spectrum_data[!, :Error]
    l = (y .- 2*y_error_array)
    u = (y .+ 2*y_error_array)
    clrs = [1, :red]
    c_confidence = clrs[2]
    α = 0.5
    plot(x, y, fillrange = u, fillalpha = α, c = c_confidence, label = "Confidence band (2σ)", seriestype=:step,
        xlims=(0,400), xticks= (0:100:500, [string(round(adc_to_energy(i), digits=0)) for i in 0:100:500]),
        xlabel="Energy (keV)", linewidth=0,
        title="MCA coincidence spectrum - Calibrated", ylabel="Detected counts",
        size=(700,500), dpi=700, titlefontsize=11, xguidefontsize=10, yguidefontsize=10,
        margin = 5mm
    )
    plot!(x, y, fillrange = l, st=:step, fillalpha = α, c = c_confidence, label = :none, linewidth=0)
    plot!(x, y, color=clrs[1], st=:step, legend=:topleft, label="Fitted line", linewidth=0.7)
    # savefig("plots\\3.4.2_MCA_coinc_calib.pdf")
end
