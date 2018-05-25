using PyPlot
include("GS.jl")
include("func.jl")

crd      = Cord(Rlen = 512, μlen = 64, a = 0.99, rmax = 100.)
mtr      = Geom(crd)
U, Ω_I, U_H   = Init(crd, mtr)
grd      = Grid(crd, mtr, Ω_I)
ils      = LS(U, grd, crd, Ω_I)
lsn      = LS_neighbors(U, ils, grd, crd)
δU       = zeros(crd.μlen)


for Ωloop = 1:10
    for Iloop = 1:5
        U, U_H, Res, dU = Solver!(U, crd, grd, Ω_I, ils, lsn, maxitr = 2, omega = 0.8)
        ils, Ω_I, δU    = IIp_updater!(U, crd, Ω_I, ils, lsn, Isf = 0.06)   #update IIp in ils and Ω_I
        grd             = Grid!(grd, crd, mtr, Ω_I)                        #update IIp in grd
        println("Iloop = $Iloop, res = $(sum(abs(Res))), U_H = $U_H")
        # plot(ils.ULS, δU)
        # plot(ils.ULS, ils.IIp/10)
    end
    Ω_I      = ΩI_updater!(U, crd, Ω_I, ils)
    grd      = Grid(crd, mtr, Ω_I)
    ils, Ω_I = LS_updater!(U, grd, crd, Ω_I, ils)     #update ils and Ω_I.Ωspl from (ils.ULS, ils.Ω)
    lsn      = LS_neighbors(U, ils, grd, crd)
    println("Ωloop = $Ωloop")
    Ubm = linspace(0., 1.1U_H, 1024)
    plot(Ubm, Ω_I.Ωspl(Ubm))
    plot(Ubm, Ω_I.Ispl(Ubm)/U_H)
    plot(Ubm, Ω_I.IIpspl(Ubm)/U_H)
end

plot(U[1, 145:154])
plot(U[2, 145:154], "--")
plot(U[3, 145:154], "k")
plot(U[4, 145:154], "k--")


Ubm = linspace(0., U_H, 2048)
fig = figure(figsize=(8,10))
subplot(311)
plot(Ubm/U_H, Ω_I.Ωspl(Ubm)/crd.Ω_H, lw = 3, "k")
ylim(0., 0.51)
ylabel(L"$\Omega/ \Omega_{\rm H}$", fontsize = 20)
tick_params(axis="both", which="major", labelsize=14)

subplot(312)
plot(Ubm/U_H, Ω_I.IIpspl(Ubm)/U_H, lw = 3, "k")
ylim(-0.025, 0.025)
ylabel(L"$II'/ A_\phi^{\rm H}$", fontsize = 20)
tick_params(axis="both", which="major", labelsize=14)

subplot(313)
plot(Ubm/U_H, 2*Ω_I.Ωspl(Ubm).*Ubm/U_H, lw = 3, "k")
ylim(0., 0.12)
xlabel(L"$A_\phi/ A_\phi^{\rm H}$", fontsize = 20)
ylabel(L"$I/ (\Omega_{\rm H}A_\phi^{\rm H})$", fontsize = 20)
tick_params(axis="both", which="major", labelsize=14)
tight_layout()
savefig("f2.pdf")



# r, fsq, favg = Fsq(U, crd, grd, lsn)
# plot(r, fsq, lw = 2)
# figure(figsize=(6,5))
# xlabel(L"$r/M$", fontsize = 20)
# ylabel(L"$\frac{B^2-E^2}{B^2+E^2}|_{\mu = 0}$", fontsize = 20)
# tick_params(axis="both", which="major", labelsize=14)
# xmajorLocator = matplotlib[:ticker][:MultipleLocator](0.2)
# ymajorLocator = matplotlib[:ticker][:MultipleLocator](0.25)
# xminorLocator = matplotlib[:ticker][:MultipleLocator](0.04)
# yminorLocator = matplotlib[:ticker][:MultipleLocator](0.05)
#
# ax = axes()
# ax[:xaxis][:set_major_locator](xmajorLocator)
# ax[:yaxis][:set_major_locator](ymajorLocator)
# ax[:xaxis][:set_minor_locator](xminorLocator)
# ax[:yaxis][:set_minor_locator](yminorLocator)
# tight_layout()
IIp_updater
