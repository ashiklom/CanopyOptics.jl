################################################################################
# using WGLMakie
f = Figure()
ax = Axis(f[1, 1])
lines!(wavelength, refl)
lines!(wavelength, 1.0 .- trans)
##########
using UnicodePlots
lineplot(wavelength, refl)
