using Dierckx
using LsqFit


function Solver!(U::Array{Float64,2}, crd::Cord, grd::Grid, Ω_I::Ω_and_I, ils::LS, lsn::LS_neighbors; maxitr = 5, omega = 0.5, ϵ = 1.0e-5)
    a = grd.aa
    b = grd.bb
    c = grd.cc
    d = grd.dd
    ee= grd.ee
    f = grd.ff

    ee_rgl = grd.ee_rgl

    Rlen     = crd.Rlen
    μlen     = crd.μlen
    idx_xbd  = crd.idx_xbd

    lsn_idx = lsn.lsn_idx
    lsn_map = lsn.lsn_map
    U_H  = 0.
    Res  = zeros(μlen, Rlen)
    dU   = zeros(U)
    Upls = zeros(crd.μlen)
    Umns = zeros(crd.μlen)

    #ρ_jcb  = ( cos(pi/μlen) * crd.δμ^2 + cos(pi/Rlen) * crd.δR^2 ) / (crd.δμ^2 + crd.δR^2)
    #ρ2_jcb = ρ_jcb.^2                                       #estimating ρ2_jcb is hard
    #ρ2_jcb = 0.05

    for n = 1:maxitr
        Res     = zeros(μlen, Rlen)
        dU      = zeros(U)

        for j = 2:μlen-1                                 #even run for even j, odd run for odd j
            for l = 2+mod(j, 2):2:idx_xbd[j]-1
                if lsn_map[j,l] == 0
                    Res[j,l] = a[j,l]*U[j+1,l] + b[j,l]*U[j-1,l] + c[j,l]*U[j,l+1] + d[j,l]*U[j,l-1] + ee[j,l]*U[j,l] - f[j,l]
                    dU[j,l]  = -omega*Res[j,l]/ee_rgl[j,l]
                end
            end
        end
        #omega = (n==1)? 1./(1.- ρ2_jcb/2.) : 1./(1.-omega*ρ2_jcb/4.)

        for j = 2:μlen-1                                 #odd run for even j, even run for odd j
            for l = 3-mod(j,2):2:idx_xbd[j]-1
                if lsn_map[j,l] == 0
                    Res[j,l] = a[j,l]*U[j+1,l] + b[j,l]*U[j-1,l] + c[j,l]*U[j,l+1] + d[j,l]*U[j,l-1] + ee[j,l]*U[j,l] - f[j,l]
                    dU[j,l]  = -omega*Res[j,l]/ee_rgl[j,l]
                end
            end
        end
        #omega = 1./(1.-omega*ρ2_jcb/4.)


        # tmp = deepcopy(dU)
        # for j = 2:crd.μlen - 1
        #     for l = 2:idx_xbd[j]-1
        #         dU[j,l] = 0.5*tmp[j,l] + 0.5*(tmp[j+1,l] + tmp[j-1, l] + tmp[j,l+1] + tmp[j,l-1])/4
        #     end
        # end

        U += dU
        U, U_H, dU  = Bounds!(U, dU, crd, Ω_I, ils, lsn)
    end
    return U, U_H, Res, dU
end

function Bounds!(U::Array{Float64,2}, dU::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I, ils::LS, lsn::LS_neighbors)
    #lsn bounds
    U = USmooth!(U, lsn, crd)                           #smooth the neighbors before interpolation

    #horizon and inf r boundary values
    U[:,1]   = U[:,2]
    U[:,end] = U[:,end-1]         #inf r boundary is not used due to xbd

    #equator boundary values (beyond ils and within)
    idx_r2  = crd.idx_r2
    idx_bd  = crd.idx_xbd[1]

    U[1, idx_r2+1:idx_bd] = U[2, idx_r2+1:idx_bd]       # ∂μ = 0, for r > 2

    Ispl = I_solver(Ω_I)
    Ωspl = Ω_I.Ωspl
    Uhe  = U[1,1]                       #U in the horizon/equator cornor
    Ωhe  = Ωspl(Uhe)
    Ihe  = Ispl(Uhe)

    Ω_H  = crd.Ω_H
    rmin = crd.rmin

    ∂μ    = zeros(idx_r2)
    ∂μ[1] = 0.5*rmin* Ihe/(Ωhe-Ω_H)      # obtain from Znajek Condition, ∂μ shoule be negative
    ∂μ[1:idx_r2] = ∂μ[1]*(crd.rcol[1:idx_r2]-2.0)/(rmin-2.0)  #initial guess: ∂μ ∈ [∂μ[1], 0] linear in r

    U[1, 1:idx_r2] = U[2, 1:idx_r2] - ∂μ[1:idx_r2]*crd.δμ       #r < 2

    δr  = crd.rcol[idx_r2+1]-crd.rcol[idx_r2]
    U_H = U[1, idx_r2]*(crd.rcol[idx_r2+1]-2.0)/δr + U[1, idx_r2+1]*(2.0-crd.rcol[idx_r2])/δr
    return U, U_H, dU
end

function Proj(U::Array{Float64,2}, crd::Cord, ils::LS, lsn::LS_neighbors)

    Rlt = lsn.ILS_lt[:,1]
    μlt = lsn.ILS_lt[:,2]
    xlt = lsn.ILS_lt[:,3]

    Rrt = lsn.ILS_rt[:,1]
    μrt = lsn.ILS_rt[:,2]
    xrt = lsn.ILS_rt[:,3]

    Uspl = Spline2D( crd.μcol, crd.Rcol, U, kx = 1, ky = 1)
    Ult  = Uspl(μlt, Rlt)
    Urt  = Uspl(μrt, Rrt)
    Umns = Ult + xlt .* ils.S
    Upls = Urt + xrt .* ils.S

    Umns[end]  = 0.
    Upls[end]  = 0.

    return Umns, Upls
end


function USmooth!(U::Array{Float64,2}, lsn::LS_neighbors, crd::Cord)
    for μidx = 2:crd.μlen-1
        Ridx = lsn.lsn_idx[μidx]
        x    = [crd.R[μidx, Ridx-2], crd.R[μidx, Ridx+3]]
        y    = [U[μidx, Ridx-2], U[μidx, Ridx+3]]        #the next near points
        spl  = Spline1D(x, y, k=1)
        U[μidx, Ridx-1] = spl(crd.R[μidx, Ridx-1])
        U[μidx, Ridx]   = spl(crd.R[μidx, Ridx])
        U[μidx, Ridx+1] = spl(crd.R[μidx, Ridx+1])
        U[μidx, Ridx+2] = spl(crd.R[μidx, Ridx+2])
    end
    return U
end

function I_updater!(U, crd, Ω_I, ils, lsn; Isf = 0.04, xbd = 4.0)
    #update IIp
    U = USmooth!(U, lsn, crd)                       #smooth the neighbors before interpolation
    Umns, Upls = Proj(U, crd, ils, lsn)

    Uils = 0.5 * (Upls + Umns)
    δU   = Upls - Umns; δU[1] = 0.                  # we expect ∂ᵣU = 0 where the LS crossing the equator

    Umodel(x, p) = (1-x) .* ( p[1]         + p[2] .* x    + p[3] .* x.^2 + p[4] .* x.^3
                            + p[5] .* x.^4 + p[6] .* x.^5 + p[7] .* x.^6 + p[8] .* x.^7)
    Ufit = curve_fit(Umodel, crd.μcol, Uils, [4., 0., 0., 0., 0., 0., 0., 0.])
    Uils = Umodel(crd.μcol, Ufit.param)

    Inew = ils.I - Isf * δU
    Imodel(x, p) =  x .* (Uils[1]-x) .* ( p[1] + p[2] .* x    + p[3] .* x.^2 + p[4] .* x.^3)
    Ifit = curve_fit(Imodel, Uils, Inew, [0., 0., 0., 0.])
    Inew = Imodel(Uils, Ifit.param)

    Ubm = collect(linspace(0, Uils[1], 1024))
    Ibm = Imodel(Ubm, Ifit.param)
    Ispl = Spline1D(Ubm, Ibm, bc = "zero")
    Ipbm = derivative(Ispl, Ubm)
    IIpspl  = Spline1D(Ubm, Ibm.*Ipbm, bc = "zero")

    #update ils and lsn (only Ispl and IIpspl related part)
    ils     = LS!(ils, Uils, Ispl, IIpspl)
    Ω_I     = Ω_and_I!(U, Ω_I, Ispl, IIpspl)
    return ils, Ω_I, δU
end




function Ω_updater!(U::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I; xbd = 4.)
    Ubm  = collect(linspace(0., xbd^2, 2048))
    Ispl = Ω_I.Ispl
    Ωbm  = Ispl(Ubm) ./ (2Ubm); Ωbm[1] = 0.5*crd.Ω_H
    Ωold = Ω_I.Ωspl(Ubm)
    Ωnew = Ωold + 0.1*(Ωbm-Ωold)
    Ωspl = Spline1D(Ubm, Ωnew)

    return Ω_and_I(U, crd, Ωspl, Ispl, Ω_I.IIpspl)
end

function Init(crd::Cord, mtr::Geom; xbd = 4.0)

        z = crd.r .* crd.μ
        x = sqrt(crd.r.^2 - z.^2)
        U = x.^2 + 1.5*acos(crd.μ).*exp(-crd.r.^2); U_H = U[1, crd.idx_r2]

        #initialize Ω_and_I
        Ubm = collect(linspace(0., xbd, 2048)).^2
        Ωbm = zeros(Ubm)
        for i = 1:length(Ubm)
            Ωbm[i] = (Ubm[i] < U_H) ? 0.5*crd.Ω_H*(cos(pi/2*Ubm[i]/U_H).^2) : 0.
        end
        Ibm = 2*Ωbm.*Ubm

        Ωspl = Spline1D(Ubm, Ωbm, k =1, bc = "zero")
        Ispl = Spline1D(Ubm, Ibm, k =1, bc = "zero")
        Ipbm  = derivative(Ispl, Ubm)
        IIpspl= Spline1D(Ubm, Ibm.*Ipbm, k =1, bc = "zero")
        Ω_I   = Ω_and_I(U, crd, Ωspl, Ispl, IIpspl)

        return U, U_H, Ω_I
end


function Rμ2xy(crd, U, ils; xmax = 3., ymax = 4., len = 512, Umax = 9.0, cnum = 30)
    spl = Spline2D( crd.μcol, crd.Rcol, U)
    rmin= 1. + sqrt(1-crd.a^2)

    x = linspace(0., xmax, len)
    y = linspace(0., ymax, len)
    x, y = np.meshgrid(x,y)
    r = sqrt(x .^2 + y .^2)
    Θ = angle(y + x .*im)

    Uxy = zeros(r)
    for i = 1:len
      for j = 1:len
        if r[i,j] > rmin
          Uxy[i,j] = evaluate(spl, cos(Θ[i,j]), r2R(r[i,j]))
        end
      end
    end

    rILS = R2r(ils.Loc[:,1])
    μILS = ils.Loc[:,2]
    yILS = rILS .* μILS
    xILS = sqrt(rILS.^2 - yILS.^2)

    rIRS = 1.+sqrt(1-crd.a^2 * μILS.^2)
    yIRS = rIRS .* μILS
    xIRS = sqrt(rIRS.^2 - yIRS.^2)

    rhz = crd.rmin*1.01
    yhz = rhz * μILS
    xhz = sqrt(rhz^2 - yhz.^2)

    levels = linspace(0.005, Umax, cnum)
    figure(figsize=(5,6))
    contour(Uxy, levels, extent = (0, xmax, 0, ymax), colors = "k")
    plot(xILS, yILS,  "k--")
    #plot(xIRS, yIRS, lw = 3, "k--")
    fill_between(xhz, 0., yhz, color = "black")
    xlabel(L"$X/M$", fontsize = 20)
    ylabel(L"$Z/M$", fontsize = 20)
    tight_layout()
    savefig("f1.pdf")
end

function Fsq(U::Array{Float64, 2}, crd::Cord, grd::Grid, lsn::LS_neighbors)
    idx = lsn.lsn_idx[1,1]
    ∂1U = (U[1, 2:idx+1] - U[1, 1:idx]) ./ crd.δR
    ∂2U = (U[2, 1:idx] - U[1, 1:idx]) ./ crd.δμ

    Rcol = crd.R[1, 1:idx]
    ∂rU  = ∂1U .* (1-Rcol).^2; ∂μU = ∂2U

    μ = 0.
    r = crd.rcol[1:idx]
    Δ = r.^2 - 2r + crd.a^2
    Σ = r.^2 + crd.a^2 * μ.^2
    β = Δ .* Σ + 2r .*(r.^2+crd.a^2)

    Ispl = I_solver(Ω_I)
    Icol = evaluate(Ispl, U[1,1:idx])
    κcol = reshape(grd.κ[1,1:idx], idx)

    B2mE2 = -κcol .* (Δ .* ∂rU.^2 + ∂μU.^2) + Σ .* Icol.^2
    B2pE2 = (κcol  + Δ .* Σ ./β ).* (Δ .* ∂rU.^2 + ∂μU.^2) +  Σ .* Icol.^2
    fsq   = B2mE2./B2pE2

    wspl = Spline1D(r[2:end], sqrt( (Σ.*β./Δ)[2:end]) )
    fspl = Spline1D(r[2:end], sqrt( (Σ.*β./Δ)[2:end]) .* fsq[2:end])
    favg = integrate(fspl, r[2], r[end]) / integrate(wspl, r[2], r[end])
    return r, fsq, favg
end

function I_solver(Ω_I::Ω_and_I; xbd = 4.0)
    Ubm = collect(linspace(0., xbd^2, 2048))
    δU  = (Ubm[end]-Ubm[1])/(length(Ubm)-1)
    IIp = Ω_I.IIpspl(Ubm)
    Isq = zeros(Ubm)

    for iter = 2:length(Isq)
        Isq[iter] = Isq[iter-1] + 2*(IIp[iter]+IIp[iter-1])*0.5*δU
    end

    Ispl = Spline1D(Ubm, sqrt(abs(Isq)), bc = "zero")
    return Ispl
end

# function cplot(ils, lsn)
#     plot(ils.Loc[:,1], ils.Loc[:,2])
#     plot(lsn.ILS_lon[:,1], lsn.ILS_lon[:,2], ".")
#     plot(lsn.ILS_ron[:,1], lsn.ILS_ron[:,2], ".")
# end
#
# function Znajek(crd::Cord, Ω_I::Ω_and_I, U_H::Float64)
#     a   = crd.a
#     Ω_H = crd.Ω_H
#     μcol= crd.μcol
#
#     Ispl = I_solver(Ω_I)
#     Ωspl = Ω_I.Ωspl
#
#     Utmp = linspace(0., U_H, 512)
#     Ispl_nm = Spline1D(Utmp/U_H, Ispl(Utmp)/U_H)
#     Ωspl_nm = Spline1D(Utmp/U_H, Ωspl(Utmp))
#
#     rmin = 1. + sqrt(1. - a^2)
#     Gμ   = (1-μcol)./(1+μcol) .* exp(a^2/rmin * μcol)
#
#     A  = collect(linspace(1., 0., 1024))
#     δA = (A[end]-A[1])/(length(A)-1)
#     IA = Ispl_nm(A); IA[end] = 0.
#     fA = IA ./ 2 ./ (Ω_H - Ωspl_nm(A))
#
#     tmp = zeros(A)
#     for i = 2:length(fA)
#         tmp[i] =  0.5*(1/fA[i-1] + 1/fA[i])
#     end
#     FA = exp(cumsum(tmp) .* δA)
#
#     Aspl = Spline1D(reverse(FA), reverse(A))
#     U_bc = Aspl(Gμ)
#     return U_H*U_bc
# end
