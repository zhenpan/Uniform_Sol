using PyPlot
using JLD
include("scripts/GS.jl")
include("scripts/SOR.jl")
include("scripts/func.jl")

Ω_sim = readdlm("/home/zhenpan/Astroph/Uniform_Sol/psi_omegaf_a_0.990.out")
@load "/home/zhenpan/Astroph/Uniform_Sol/scripts/a99_Icrt.jld"

Ubm = linspace(0., U_H, 64)
fig = figure(figsize=(5, 3.6))
plot(Ubm/U_H, Ω_I.Ωspl(Ubm)/crd.Ω_H, lw = 2, "k")
plot(Ubm/U_H, 0.5*(1-Ubm/U_H), lw = 2, "r--")
plot(Ω_sim[:,1]/Ω_sim[end,1], Ω_sim[:,2], lw = 2)
ylim(0., 0.52)
xlabel(L"$A_\phi$", fontsize = 20)
ylabel(L"$\Omega/ \Omega_{\rm H}$", fontsize = 20)
tick_params(axis="both", which="major", length =4, labelsize=14)
ax = gca()
mx = matplotlib[:ticker][:MultipleLocator](0.2)
My = matplotlib[:ticker][:MultipleLocator](0.1) # Define interval of major ticks
my = matplotlib[:ticker][:MultipleLocator](0.05)
ax[:xaxis][:set_minor_locator](mx)
ax[:yaxis][:set_major_locator](My)
ax[:yaxis][:set_minor_locator](my)
tight_layout()
savefig("f2b.pdf")
