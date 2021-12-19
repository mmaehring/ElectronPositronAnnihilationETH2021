# using Plots
begin
   using StatsPlots
   using Measurements
   using Unitful
   using UnitfulRecipes
end

voltages = [0, 50, 100, 150, 200, 300, 400, 500, 600, 625, 650, 675, 700,
           725,  750, 760,  775, 780, 800,  825,  850,  875, 900, 910,
           925,  950,  975, 1000, 1050, 1100, 1150, 1200] # V

counts = [0.0, 0.0 , 1.0  , 0.0 ,  0.0,   0.0,   0.0,   2.0,  67.0, 187.0, 321.0, 879.0, 1218.0, 1664.0,
          1836.0, 1909.0, 1920.0, 1950.0, 2064.0, 2210.0, 2206.0, 2239.0, 2336.0, 2373.0, 2412.0, 2493.0,
          2783.0, 3154.0, 7600.0, 86793.0, 256382.0, 798189.0] # counts

counts_errors = 2*sqrt.(counts) # 2 σ

count_and_errors = counts .± counts_errors

begin
   plot(voltages[8:end], counts[8:end], yaxis=:log,
   xticks=[500, 550, 600, 650, 700, 750, 800,  850, 900,  950, 1000, 1050, 1100, 1150, 1200],
   xlabel="HV (V)", ylabel="Measured counts (log)",
   title="Working point inspection of the PMT",
   label=:none, legend=:bottomright,
   titlefontsize = 8
   );
   scatter!(voltages[9:end], count_and_errors[9:end], label="Measured points ± 2σ", yaxis=:log, markersize=2, ylims=(1e-0, 1e6));
   scatter!([voltages[8]], [counts[8] ± 0.8], label=:none, yaxis=:log, color=:orange, markersize=2, ylims=(1e-0, 1e6));
   vline!([850], label="Suggested working point", color=:black);
   plot!(500:100:1250, 2162*ones(length(500:100:1250)),
      xlims=(485, 1225), color=:cyan,
      linestyle=:dash, # ribbon=(221.65),
      label="Counts at suggested working point", legendfont=font(7)
   );
   annotate!(720,5000, text("Average count of the values in \n the suggested range: 2162 ± 220", 7, :dark))
end

# savefig("plots\\working_point_of_pmt.png")
