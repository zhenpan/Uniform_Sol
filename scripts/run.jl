using PyPlot
using JLD
include("scripts/GS.jl")
include("scripts/SOR.jl")
include("scripts/func.jl")

crd      = Cord(Rlen = 512, μlen = 64, a = 0.99, rmax = 100., xbd = 4.5)
mtr      = Geom(crd)
Ω_par    = [0.1, 0.02, 0.04, 0.0]
U, Ω_I, U_H   = Init(crd, mtr, Ω_par)
grd      = Grid(crd, mtr, Ω_I)
ils      = LS(U, grd, crd, Ω_I)
lsn      = LS_neighbors(U, ils, grd, crd)
bc_eqt   = BC_gen(U, crd, grd, Ω_I, BC_opt = 0)  # 0:init, else:updater
Ωpar_loop= 1

δU = zeros(crd.μlen); Res = 0.; fsq2_avg  = 0.1

while (Ωpar_loop <= 15) & (fsq2_avg > 2.e-4)
    for BCloop = 1:15
        for Ωloop = 1:150
            for Iloop = 1:5
                U,  Res, dU = Solver!(U, crd, grd, Ω_I, ils, lsn, maxitr = 2, omega = 1.0)
                U, U_H      = Bounds!(U, crd, Ω_I, lsn, bc_eqt)                  #where bc_eqt contains ∂μU on equator
                ils,Ω_I, δU = IIp_updater!(U, crd, Ω_I, ils, lsn, Isf = 0.08)    #update IIp in ils and Ω_I
                grd         = Grid!(grd, crd, mtr, Ω_I)                          #update IIp in grd
            end
            println("(Ωpar-BC-Ωloop) = ($Ωpar_loop $BCloop $Ωloop), Res = $(sum(abs(Res))), U_H = $U_H, U_he = $(U[1,1])")
            Ω_I      = ΩI_updater!(U, crd, Ω_I, ils, Ω_par)
            grd      = Grid(crd, mtr, Ω_I)
            ils, Ω_I = LS_updater!(U, grd, crd, Ω_I, ils, Ω_par)     #update ils and Ω_I.Ωspl from (ils.ULS, ils.Ω)
            lsn      = LS_neighbors(U, ils, grd, crd)
        end
            bc_eqt   = BC_gen(U, crd, grd, Ω_I, BC_opt = 1, Isf = 6.0)
            Ubm = linspace(0., U_H, 64)
            plot(Ubm, Ω_I.Ωspl(Ubm))
            plot(Ubm, Ω_I.Ispl(Ubm)/U_H)
            plot(Ubm, Ω_I.IIpspl(Ubm)/U_H)

            Ispl = Spline1D(Ubm, 2*Ω_I.Ωspl(Ubm).*Ubm)
            plot(Ubm, 0.5*crd.Ω_H*(1-Ubm/U_H), "k--")
            plot(Ubm, Ispl(Ubm)/U_H, "k--")
            plot(Ubm, Ispl(Ubm).*derivative(Ispl, collect(Ubm))/U_H, "k--")
    end
    Ω_par, fsq2_avg =  Ωpar_updater!(U, crd, grd, Ω_I, ils, Isf = 0.06)
    Ωpar_loop +=1
    println("fsq2_avg = $fsq2_avg")
    @save "/home/zhenpan/Astroph/Uniform_Sol/scripts/a$(crd.a)_test.jld" Ωpar_loop crd mtr Ω_par U Ω_I U_H grd ils lsn δU bc_eqt fsq2_avg
end


# plot(U[1, 145:154])
# plot(U[2, 145:154], "--")
# plot(U[3, 145:154], "k")
# plot(U[4, 145:154], "k--")

Ueqt, B2mE2, B2mE2_expt, fsq, fsq2_avg = Fsq(U, crd, grd, Ω_I)
plot(Ueqt, fsq, lw = 2)
plot(Ueqt, zeros(Ueqt), "--")

######################################################################
#reload
######################################################################


using JLD
using PyPlot
include("scripts/GS.jl")
include("scripts/SOR.jl")
include("scripts/func.jl")

@load "/home/zhenpan/Astroph/Uniform_Sol/scripts/a0.99.jld"

δU = zeros(crd.μlen); Res = 0.

while (Ωpar_loop <= 45) & (fsq2_avg > 1.e-4)
    for BCloop = 1:15
        for Ωloop = 1:150
            for Iloop = 1:5
                U,  Res, dU = Solver!(U, crd, grd, Ω_I, ils, lsn, maxitr = 2, omega = 1.0)
                U, U_H      = Bounds!(U, crd, Ω_I, lsn, bc_eqt)                  #where bc_eqt contains ∂μU on equator
                ils,Ω_I, δU = IIp_updater!(U, crd, Ω_I, ils, lsn, Isf = 0.08)    #update IIp in ils and Ω_I
                grd         = Grid!(grd, crd, mtr, Ω_I)                          #update IIp in grd
            end
            println("(Ωpar-BC-Ωloop) = ($Ωpar_loop $BCloop $Ωloop), Res = $(sum(abs(Res))), U_H = $U_H, U_he = $(U[1,1])")
            Ω_I      = ΩI_updater!(U, crd, Ω_I, ils, Ω_par)
            grd      = Grid(crd, mtr, Ω_I)
            ils, Ω_I = LS_updater!(U, grd, crd, Ω_I, ils, Ω_par)     #update ils and Ω_I.Ωspl from (ils.ULS, ils.Ω)
            lsn      = LS_neighbors(U, ils, grd, crd)
        end
            ISF_BE   = Ωpar_loop > 10 ? 0.2:0.
            bc_eqt   = BC_gen(U, crd, grd, Ω_I, BC_opt = 1, Isf = 6.0, Isf_BE = ISF_BE)
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
    @save "/home/zhenpan/Astroph/Uniform_Sol/scripts/a$(crd.a).jld" Ωpar_loop crd mtr Ω_par U Ω_I U_H grd ils lsn δU bc_eqt fsq2_avg
end



###############################################################################

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
