using PyPlot
include("GS.jl")
include("func.jl")

crd      = Cord(Rlen = 512, μlen = 64, a = 0.99, rmax = 100.)
mtr      = Geom(crd)
U_H      = 4.5
U, Ω_I   = Init(crd, mtr, U_H = U_H)
grd      = Grid(crd, mtr, Ω_I)
ils      = LS(U, grd, crd, Ω_I)
lsn      = LS_neighbors(U, ils, grd, crd)


for Ωloop = 1:200
    for Iloop = 1:10
        U, U_H, Res, dU = Solver!(U, crd, grd, ils, lsn, maxitr = 2, omega = 0.8)
        ils, Ω_I, δU    = I_updater!(U, crd, Ω_I, ils, lsn, Isf = 0.02)   #update IIp in ils and Ω_I
        grd             = Grid!(grd, crd, mtr, Ω_I)           #update IIp in grd
        println("Iloop = $Iloop, res = $(sum(abs(Res))), U_H = $U_H")
    end
    # Ω_I   = Ω_updater!(U, U_H, crd, Ω_I)
    Ω_I   = Ω_updater!(U, crd, Ω_I)
    grd   = Grid(crd, mtr, Ω_I)
    ils   = LS(U, grd, crd, Ω_I)
    lsn   = LS_neighbors(U, ils, grd, crd)
    println("loop = $Ωloop")
    plot(1.1ils.ULS/U_H, Ω_I.IIpspl(1.1ils.ULS)/U_H)
    plot(1.1ils.ULS/U_H, Ω_I.Ωspl(1.1ils.ULS))
end

Ubm = linspace(0., U_H, 256)
subplot(211)
plot(Ubm/U_H, Ω_I.Ωspl(Ubm)/crd.Ω_H, lw = 3, "k")
ylim(0.28, 0.51)
ylabel(L"$\Omega/ \Omega_{\rm H}$", fontsize = 20)
tick_params(axis="both", which="major", labelsize=14)

subplot(212)
plot(Ubm/U_H, Ω_I.IIpspl(Ubm), lw = 3, "k")
ylim(0., 0.21)
xlabel(L"$A_\phi/ A_\phi^{\rm H}$", fontsize = 20)
ylabel(L"$II'/ A_\phi^{\rm H}$", fontsize = 20)
tick_params(axis="both", which="major", labelsize=14)
tight_layout()
savefig("f2.pdf")




#plot(Ubm/U_H, sqrt(1-Ubm/U_H/1.1)./(1.+sqrt(1-Ubm/U_H/1.1)))

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
