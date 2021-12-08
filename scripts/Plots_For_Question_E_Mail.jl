using Plots
using StatsPlots

voltages = [0, 50, 100, 150, 200, 300, 400, 500, 600, 625, 650, 675, 700,
           725,  750, 760,  775, 780, 800,  825,  850,  875, 900, 910,
           925,  950,  975, 1000, 1050, 1100, 1150, 1200] # V

counts = [0, 0 , 1  , 0  ,   0,   0,   0,   2,  67, 187, 321, 879, 1218, 1664,
          1836, 1909, 1920, 1950, 2064, 2210, 2206, 2239, 2336, 2373, 2412, 2493,
          2783, 3154, 7600, 86793, 256382, 798189] # counts

plot(voltages[8:end], counts[8:end], yaxis=:log,
   xticks=[500, 550, 600, 650, 700, 750, 800,  850, 900,  950, 1000, 1050, 1100, 1150, 1200],
   xlabel="HV (V)", ylabel="Amount of measured counts (log)",
   title="Working point inspection of the PMT",
   label=:none, legend=:bottomright,
   titlefontsize = 8
);
scatter!(voltages[8:end], counts[8:end], label="Measured points");
vline!([850], label="Suggested working point");
plot!(500:100:1250, 2162*ones(length(500:100:1250)),
   xlims=(485, 1225),
   linestyle=:dash, # ribbon=(221.65),
   label="Counts at suggested working point", legendfont=font(7)
);
annotate!(720,5000,
         Plots.text("Average count of the values in \n the suggested range: 2162 ± 220", 7, :dark))
# savefig("working_point_of_pmt.png")
