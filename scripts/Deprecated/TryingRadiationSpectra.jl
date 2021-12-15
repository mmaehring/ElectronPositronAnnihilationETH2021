## Trying out RadiationSpectra
# Finding peaks & calibration
# Utility for radiation spectra
using RadiationSpectra
h_uncal = RadiationSpectra.get_example_spectrum()
photon_lines = [609.312, 911.204, 1120.287, 1460.830, 1764.494]; # keV
h_cal, h_deconv, peakPositions, threshold, c, c_precal = RadiationSpectra.calibrate_spectrum(h_uncal, photon_lines)

p_uncal = plot(h_uncal, st=:step, label="Uncalibrated spectrum");
p_deconv = plot(h_deconv, st=:step, label = "Deconvoluted spectrum");
hline!([threshold], label = "threshold", lw = 1.5);
p_cal = plot(h_cal, st=:step, label="Calibrated spectrum", xlabel="E / keV");
vline!(photon_lines, lw=0.5, color=:red, label="Photon lines");
plot(p_uncal, p_deconv, p_cal, size=(800,700), layout=(3, 1));
