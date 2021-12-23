## Advanced fit -> Leo 1994  7.2 -> pg 158
# N_7_3(t, p) = @. p[1] * exp(-t/p[2]) + p[3] * exp(-t/p[4])
#
# p0 = [500.0, 230.0, 500.0, 230.0];
# lb = zeros(length(p0));
#
# fit2 = curve_fit(N_7_3, ustrip.(t), ustrip.(h), ustrip(w), p0, lower=lb);
# adv_params = fit2.param
# adv_errrors = stderror(fit2) # LAPACK exception
#
# fitted_function(t) = N_7_3(ustrip.(t), params)
#
# p_adv = begin
#     plot(t .+ df[!, :Time][strt], h, color=:blue, label=:none);
#     plot!(df[!, :Time][begin:strt], df[!, :Height][begin:strt], color=:blue, label="Collected data"); # plot data
#
#     annotate!(375, 15, text(L"N = A \cdot \exp\left(\frac{-t}{\tau_f}\right) \! + \! B \cdot \exp\left(\frac{-t}{\tau_s}\right)", :red, 10))
#     plot!(t .+ df[!, :Time][strt], N_7_3(ustrip.(t), ustrip.(adv_params)),
#             color=:red, label= "Fit to data, decay time: $((adv_params[2] ± 8.9)*u"ns")",
#             xlabel="Time series",
#             ylabel="Height (mV)",
#             legendfontsize=7,
#             axisfontsize=4
#     ) # add fit
#     # scatter!(df[!, :Time][1:every:end], (df[!, :Height] .± df[!, :σ])[1:every:end], label="Selected errors", markersize=3)
#     # savefig(save_loc)
# end

## plot both parts
#
# l = @layout [a{0.5h};
#              b{0.5h}]

# plot(p_ez, p_adv, layout=l) # dpi=900)
