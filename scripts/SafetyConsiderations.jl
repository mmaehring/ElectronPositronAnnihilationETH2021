## Activate envirnoment
# using DrWatson
# quickactive("C:\\Users\\marcu\\OneDrive\\Desktop\\PraktikumIII\\e+e-_Annihilation")
## Loading libraries
using Measurements
using Unitful
using Plots
using UnitfulRecipes

## Calculating the activity of the sample
T1_2 = 2.6
λ = log(2) / T1_2
N₂₀₀₉ = 370u"kBq"

N₂₀₂₁ = N₂₀₀₉ * exp(-12.5*log(2)/2.6)

## Plotting sample activity
t = 0.0001:0.05:15.5
T1_2 = 2.6
λ = log(2) / T1_2

exp_decay(t) = N₂₀₀₉ * exp(-λ*t)

plot(t, exp_decay.(t),
    title="Decay of the ²²Na source",
    xlims=(0,15),
    xlabel="Time (years)",
    ylabel="Activity",
    label="Activity of the ²²Na sample"
)

# savefig("SampleDecay.png")

#### How much radiation penetrates the shielding?
## Beer lambert

# Got μ/ρ from: https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html
ρ = (11.34 ± 0.005)u"g/cm^3"
μ = ρ * (16.140 ± 0.005)*1e-2u"cm^2/g"
I₀ = N₂₀₂₁
I(x) = I₀ * exp(-x*μ)

intensity_after_beam = I(5u"cm") |> u"Bq"

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

5
