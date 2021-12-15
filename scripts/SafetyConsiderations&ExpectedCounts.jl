## Activate envirnoment
# using DrWatson
# quickactive("C:\\Users\\marcu\\OneDrive\\Desktop\\PraktikumIII\\e+e-_Annihilation")
## Loading libraries
using Measurements
using Unitful
using StatsPlots
using UnitfulRecipes
using Interpolations
using CSV
using DataFrames

## Calculating the activity of the sample
T1_2 = 2.601 ± 0.001
λ = log(2) / T1_2
N₂₀₀₉ = (370±0.05)u"kBq"
t = 12.5±0.05

N₂₀₂₁ = N₂₀₀₉ * exp(-t*log(2)/2.6)


## Plotting sample activity
t_linsp = 0.0001:0.05:15.5
exp_decay(t) = N₂₀₀₉ * exp(-λ*t)

plot(t_linsp, Measurements.value.(exp_decay.(t_linsp)),
    title="Decay of the ²²Na source",
    xlims=(0,15),
    xlabel="Time (years)",
    ylabel="Activity",
    label="Activity of the ²²Na sample",
)
# scatter!(t_linsp[1:10:end], exp_decay.(t_linsp[1:10:end]), label=:none, markersize=0.5)

# savefig("SampleDecay.png")

#### How much radiation penetrates the shielding?
## Beer lambert

# Got μ/ρ from: https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html -> also txt
μ_ρ_data = CSV.read("scripts\\AttenuationDataLead.txt", DataFrame)
μ_ρ_data[!,:EnergyMeV] = convert.(Float64, μ_ρ_data[!,:EnergyMeV])

ρ = (11.34 ± 0.005)u"g/cm^3"
μ511keV = ρ * (16.140 ± 0.005)*1e-2u"cm^2/g"
I₀ = N₂₀₂₁
I(x, μ) = I₀ * exp(-x*μ)

intensity_after_beam = I(5u"cm", μ511keV) |> u"Bq"

#### Determining yearly dose
h₁₀ = 0.33u"(mSv / hr) / (GBq)"
h₁₀_Bq = 0.33u"(mSv / hr) / (GBq)" |>u"(mSv / hr) / (Bq)"

dose_per_hour_unshielded = h₁₀_Bq * N₂₀₂₁ |>u"(mSv / hr)"
dose_per_hour_shielded = h₁₀_Bq * intensity_after_beam |>u"(mSv / hr)"

approx_distance = 2.5u"m"
approx_time = 12u"hr"

total_dose(dose_per_hour, distance, time) = dose_per_hour * time / ustrip(distance)^2

approx_total_dose_upper_bound = total_dose(dose_per_hour_unshielded, approx_distance, approx_time) |>u"nSv"
approx_total_dose_upper_est = total_dose(Measurements.value(dose_per_hour_shielded),
                                    approx_distance, approx_time) |>u"pSv"


## Expected counts
d = 2.54u"cm"
r = 5u"cm"
# N̄1 = (π * (d/2)^2 ) / (4π*r^2) * N₂₀₂₁ |> u"Bq"
N̄ = d^2 / (16 * r^2) * N₂₀₂₁  |> u"Bq"
# numerical_error = (abs(N̄ - N̄1), N̄/N̄1)

total_counts_upper_bound = N̄ * 10 # second measurment interval
# assume travel distance d from above
μ_NaTi_511keV = 0.3u"cm^-1"

real_detection = total_counts_upper_bound* exp(-μ_NaTi_511keV * d)
