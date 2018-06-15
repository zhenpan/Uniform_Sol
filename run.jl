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

for Ωpar_loop = 1:3
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
    Ω_par =  Ωpar_updater!(U, crd, grd, Ω_I, ils, Isf = 0.1)
end


plot(U[1, 145:154])
plot(U[2, 145:154], "--")
plot(U[3, 145:154], "k")
plot(U[4, 145:154], "k--")


Ueqt, B2mE2, fsq, fsq2_avg = Fsq(U, crd, grd, Ω_I)
plot(Ueqt, fsq, lw = 2)
plot(Ueqt, zeros(Ucol), "--")

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
