using PyPlot
include("GS.jl")
include("SOR.jl")
include("func.jl")

crd      = Cord(Rlen = 512, μlen = 64, a = 0.99, rmax = 100., xbd = 3.0)
mtr      = Geom(crd)
Ω_par    = [0.1, 0.02, 0.04, 0.0]
U, Ω_I, U_H   = Init(crd, mtr, Ω_par)
grd      = Grid(crd, mtr, Ω_I)
ils      = LS(U, grd, crd, Ω_I)
lsn      = LS_neighbors(U, ils, grd, crd)
δU       = zeros(crd.μlen); Res = 0.
bc_eqt   = BC_gen(U, crd, Ω_I, BC_opt = 0)  # 0:init, else:updater

Ωpar_loop = 1; fsq2_avg  = 0.1

while (Ωpar_loop <= 15) & (fsq2_avg > 1.e-4)
    for BCloop = 1:20
        for Ωloop = 1:200
            for Iloop = 1:5
                U,  Res, dU = Solver!(U, crd, grd, Ω_I, ils, lsn, maxitr = 2, omega = 1.0)
                U, U_H      = Bounds!(U, crd, Ω_I, lsn, bc_eqt)                  #where bc_eqt contains ∂μU on equator
                ils,Ω_I, δU = IIp_updater!(U, crd, Ω_I, ils, lsn, Isf = 0.08)    #update IIp in ils and Ω_I
                grd         = Grid!(grd, crd, mtr, Ω_I)                          #update IIp in grd
                #println("(BCloop, Ωloop, Iloop) = ($BCloop $Ωloop $Iloop), res = $(sum(abs(Res))), U_H = $U_H, U_he = $(U[1,1])")
            end
            println("(Ωpar-BC-Ωloop) = ($Ωpar_loop $BCloop $Ωloop), Res = $(sum(abs(Res))), U_H = $U_H, U_he = $(U[1,1])")
            Ω_I      = ΩI_updater!(U, crd, Ω_I, ils, Ω_par)
            grd      = Grid(crd, mtr, Ω_I)
            ils, Ω_I = LS_updater!(U, grd, crd, Ω_I, ils, Ω_par)     #update ils and Ω_I.Ωspl from (ils.ULS, ils.Ω)
            lsn      = LS_neighbors(U, ils, grd, crd)
        end
            bc_eqt   = BC_gen(U, crd, Ω_I, BC_opt = 1, Isf = 6.0)
            Ubm = linspace(0., U_H, 64)
            plot(Ubm, Ω_I.Ωspl(Ubm))
            plot(Ubm, Ω_I.Ispl(Ubm)/U_H)
            plot(Ubm, Ω_I.IIpspl(Ubm)/U_H)

            Ispl = Spline1D(Ubm, 2*Ω_I.Ωspl(Ubm).*Ubm)
            plot(Ubm, 0.5*crd.Ω_H*(1-Ubm/U_H), "k--")
            plot(Ubm, Ispl(Ubm)/U_H, "k--")
            plot(Ubm, Ispl(Ubm).*derivative(Ispl, collect(Ubm))/U_H, "k--")
    end
    Ω_par, fsq2_avg =  Ωpar_updater!(U, crd, grd, Ω_I, ils, Isf = 0.1)
    Ωpar_loop +=1
    println("fsq2_avg = $fsq2_avg")
end


# plot(U[1, 145:154])
# plot(U[2, 145:154], "--")
# plot(U[3, 145:154], "k")
# plot(U[4, 145:154], "k--")


Ueqt, B2mE2, fsq, fsq2_avg = Fsq(U, crd, grd, Ω_I)
plot(Ueqt, fsq, lw = 2)
plot(Ueqt, zeros(Ueqt), "--")

using JLD
@save "/home/zhenpan/Astroph/Uniform_Sol/a99.jld" crd mtr Ω_par U Ω_I U_H grd ils Ueqt fsq fsq2_avg
#@save "/home/zhenpan/Astroph/Uniform_Sol/a998.jld" crd mtr Ω_par U Ω_I U_H grd ils Ueqt fsq fsq2_avg

using JLD
using PyPlot
include("GS.jl")
include("SOR.jl")
include("func.jl")

@load "/home/zhenpan/Astroph/Uniform_Sol/a99.jld"

Ubm = linspace(0., U_H, 64)
fig = figure(figsize=(5, 3.6))
plot(Ubm, Ω_I.Ωspl(Ubm)/crd.Ω_H, lw = 2, "k")
plot(Ubm, 0.5*(1-Ubm/U_H), lw = 2, "r--")
ylim(0., 0.52)
xlabel(L"$A_\phi$", fontsize = 20)
ylabel(L"$\Omega/ \Omega_{\rm H}$", fontsize = 20)
tick_params(axis="both", which="major", length =4, labelsize=14)
ax = gca()
mx = matplotlib[:ticker][:MultipleLocator](0.5)
My = matplotlib[:ticker][:MultipleLocator](0.1) # Define interval of major ticks
my = matplotlib[:ticker][:MultipleLocator](0.05)
ax[:xaxis][:set_minor_locator](mx)
ax[:yaxis][:set_major_locator](My)
ax[:yaxis][:set_minor_locator](my)
tight_layout()
savefig("f2.pdf")

# fig = figure(figsize=(8,10))
# subplot(311)
# plot(Ubm/U_H, Ω_I.Ωspl(Ubm)/crd.Ω_H, lw = 2, "k")
# plot(Ubm/U_H, 0.5*(1-Ubm/U_H), lw = 2, "r--")
# ylim(0., 0.52)
# ylabel(L"$\Omega/ \Omega_{\rm H}$", fontsize = 20)
# tick_params(axis="both", which="major", labelsize=14)
# ax = gca()
# My = matplotlib[:ticker][:MultipleLocator](0.1) # Define interval of major ticks
# ax[:yaxis][:set_major_locator](My)
#
# subplot(312)
# plot(Ubm/U_H, Ω_I.IIpspl(Ubm)/(crd.Ω_H*U_H), lw = 2, "k")
# ylim(-0.07, 0.062)
# ylabel(L"$II'/ (\Omega_{\rm H}A_\phi^{\rm H})$", fontsize = 20)
# tick_params(axis="both", which="major", labelsize=14)
# bx = gca()
# My = matplotlib[:ticker][:MultipleLocator](0.03) # Define interval of major ticks
# bx[:yaxis][:set_major_locator](My)
#
# subplot(313)
# plot(Ubm/U_H, 2*Ω_I.Ωspl(Ubm).*Ubm/(crd.Ω_H*U_H), lw = 2, "k")
# ylim(0., 0.32)
# xlabel(L"$A_\phi/ A_\phi^{\rm H}$", fontsize = 20)
# ylabel(L"$I/ (\Omega_{\rm H}A_\phi^{\rm H})$", fontsize = 20)
# tick_params(axis="both", which="major",labelsize=14)
# cx = gca()
# My = matplotlib[:ticker][:MultipleLocator](0.06) # Define interval of major ticks
# cx[:yaxis][:set_major_locator](My)
#
# tight_layout()
# savefig("f2.pdf")
