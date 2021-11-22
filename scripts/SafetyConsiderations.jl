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
I(x) = 
