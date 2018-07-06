using Dierckx
using LsqFit

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
        if Ridx > 2
            x    = [crd.R[μidx, Ridx-2], crd.R[μidx, Ridx+3]]
            y    = [U[μidx, Ridx-2], U[μidx, Ridx+3]]        #the next near points
            spl  = Spline1D(x, y, k=1)
            U[μidx, Ridx-1] = spl(crd.R[μidx, Ridx-1])
            U[μidx, Ridx]   = spl(crd.R[μidx, Ridx])
            U[μidx, Ridx+1] = spl(crd.R[μidx, Ridx+1])
            U[μidx, Ridx+2] = spl(crd.R[μidx, Ridx+2])
        end
    end

    for μidx = 2:crd.μlen-1
        Ridx = lsn.lsn_idx[μidx]

        U[μidx, Ridx-1] = (U[μidx, Ridx-1] < U[μidx-1, Ridx-1]) ? U[μidx, Ridx-1]: 0.5*(U[μidx+1, Ridx-1]+U[μidx-1, Ridx-1])
        U[μidx, Ridx]   = (U[μidx, Ridx  ] < U[μidx-1, Ridx  ]) ? U[μidx, Ridx  ]: 0.5*(U[μidx+1, Ridx  ]+U[μidx-1, Ridx  ])
        U[μidx, Ridx+1] = (U[μidx, Ridx+1] < U[μidx-1, Ridx+1]) ? U[μidx, Ridx+1]: 0.5*(U[μidx+1, Ridx+1]+U[μidx-1, Ridx+1])
        U[μidx, Ridx+2] = (U[μidx, Ridx+2] < U[μidx-1, Ridx+2]) ? U[μidx, Ridx+2]: 0.5*(U[μidx+1, Ridx+2]+U[μidx-1, Ridx+2])
    end

    for μidx = 2:crd.μlen-1
        Ridx = lsn.lsn_idx[μidx]

        U[μidx, Ridx-1] = (U[μidx, Ridx-1] > U[μidx+1, Ridx-1]) ? U[μidx, Ridx-1]: 0.5*(U[μidx+1, Ridx-1]+U[μidx-1, Ridx-1])
        U[μidx, Ridx]   = (U[μidx, Ridx  ] > U[μidx+1, Ridx  ]) ? U[μidx, Ridx  ]: 0.5*(U[μidx+1, Ridx  ]+U[μidx-1, Ridx  ])
        U[μidx, Ridx+1] = (U[μidx, Ridx+1] > U[μidx+1, Ridx+1]) ? U[μidx, Ridx+1]: 0.5*(U[μidx+1, Ridx+1]+U[μidx-1, Ridx+1])
        U[μidx, Ridx+2] = (U[μidx, Ridx+2] > U[μidx+1, Ridx+2]) ? U[μidx, Ridx+2]: 0.5*(U[μidx+1, Ridx+2]+U[μidx-1, Ridx+2])
    end

    return U
end


#update IIp
function IIp_updater!(U, crd, Ω_I, ils, lsn; Isf = 0.02, xbd = 4.0)
    U = USmooth!(U, lsn, crd)                       #smooth the neighbors before interpolation
    Umns, Upls = Proj(U, crd, ils, lsn)

    Uils = 0.5 * (Upls + Umns)

    Ridx = lsn.lsn_idx[:,1]
    RILS = ils.Loc[:,1]
    Uils = zeros(RILS)
    for μidx = 1: length(RILS)
        ridx = Ridx[μidx]
        δR   = crd.δR
        Uils[μidx] = U[μidx, ridx]*(crd.R[μidx, ridx+1]-RILS[μidx])/δR + U[μidx, ridx+1]*(RILS[μidx]-crd.R[μidx, ridx])/δR
    end

    δU   = Upls - Umns; δU[1] = 0.

    IIpnew = ils.IIp -  δU*Isf; IIpnew[1] = 0.
    IIpmodel(x, p) = crd.Ω_H^2 * x .* (1.          + p[1] .* x    + p[2] .* x.^2
                                    + p[3] .* x.^3 + p[4] .* x.^4 + p[5] .* x.^5
                                    + p[6] .* x.^6 + p[7] .* x.^7 + p[8] .* x.^8)

    IIpfit = curve_fit(IIpmodel, Uils, IIpnew, [0., 0., 0., 0., 0., 0., 0., 0.])
    IIpnew = IIpmodel(Uils, IIpfit.param)

    IIpspl = IIp_gen(Uils, IIpnew)

    #update ils and lsn (only IIpspl related part)
    ils     = LS!(ils, Uils, Ω_I.Ispl, IIpspl)
    Ω_I     = Ω_and_I!(U, Ω_I, IIpspl)
    return ils, Ω_I, δU
end

function IIp_gen(Uils::Array{Float64,1}, IIp::Array{Float64,1}; drc = 0.1, xbd = 4.)

    IIpspl = Spline1D(reverse(Uils), reverse(IIp))
    U_H    = Uils[1]
    Ubm    = vcat(linspace(0., U_H, 1024), linspace(U_H*1.0001, xbd^2, 1024))
    IIpbm  = zeros(Ubm)

    for iter in eachindex(Ubm)
        if Ubm[iter] <= Uils[1]
            IIpbm[iter] = IIpspl(Ubm[iter])
        else
            IIpbm[iter] = 0.
        end
    end

    IIpspl = Spline1D(Ubm, IIpbm, k=1, bc="zero")
    return IIpspl
end

function Ω_fnc(Ω_H::Float64, Ω_par::Array{Float64}, xcol::Array{Float64})
        return Ω_H.*(1-xcol).*(0.5 + Ω_par[1]*xcol + Ω_par[2]*xcol.^2 + Ω_par[3]*xcol.^3 + Ω_par[4]*xcol.^4)     #xcol = Ucol/U_H
end

function Ωpar_updater!(U::Array{Float64,2}, crd::Cord, grd::Grid, Ω_I::Ω_and_I, ils::LS; Isf = 0.02)
    Ucol = reverse(ils.ULS); U_H =  Ucol[end]

    Ueqt, δ_B2mE2, fsq, fsq2_avg = Fsq(U, crd, grd, Ω_I)
    Icrt_spl  = Spline1D(Ueqt[2:end], δ_B2mE2[2:end], bc="zero")
    I_eqt     = 2*Ω_I.Ωspl(Ucol).*Ucol - Isf*Icrt_spl(Ucol)
    I_hrz     = Ω_I.Ispl(Ucol)
    I_new     = I_hrz.*(Ucol[end] -Ucol)/(Ucol[end]-Ucol[1]) + I_eqt.*(Ucol-Ucol[1])/(Ucol[end]-Ucol[1])

    xcol = Ucol/U_H;  Ω_H  = crd.Ω_H
    Ωmodel(x, p) = Ω_H.*(1-x).*(0.5+p[1]*x + p[2]*x.^2 + p[3]*x.^3 + p[4]*x.^4)
    Imodel(x, p) = 2*(U_H*x).*(  Ω_H.*(1-x).*(0.5+p[1]*x + p[2]*x.^2 + p[3]*x.^3 + p[4]*x.^4) )
    Ifit         = curve_fit(Imodel, xcol, I_new, [0., 0., 0., 0.])
    Ωnew         = Ωmodel(xcol, Ifit.param)
    return Ifit.param, fsq2_avg
end



function ΩI_updater!(U::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I, ils::LS, Ω_par::Array{Float64})
    U_H  = ils.ULS[1]
    Ispl = I_solver(Ω_I, U_H)
    Ωnew = Ω_fnc(crd.Ω_H, Ω_par, ils.ULS/U_H)
    Ωspl = Spline1D(reverse(ils.ULS), reverse(Ωnew), bc="zero")
    return Ω_and_I(U, crd, Ωspl, Ispl, Ω_I.IIpspl)
end

function I_solver(Ω_I::Ω_and_I, U_H::Float64)
    Ubm = linspace(0., U_H, 1024)
    δU  = U_H/(length(Ubm)-1)
    IIp = Ω_I.IIpspl(Ubm)
    Isq = zeros(Ubm)

    iter = 2
    while iter <= length(Isq) && Isq[iter-1] >= 0.
        Isq[iter] = Isq[iter-1] + 2*(IIp[iter]+IIp[iter-1])*0.5*δU
        iter = iter + 1
    end

    Isq  = max(Isq, 0.)
    Ispl = Spline1D(Ubm, sqrt(Isq), bc = "zero")

    return Ispl
end


function Rμ2xy(crd, U, ils; xmax = 3., ymax = 4., len = 1024, Umax = 9.0, cnum = 30)
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
       else
          Uxy[i,j] = evaluate(spl, cos(Θ[i,j]), r2R(rmin))
        end
      end
    end

    μcol = linspace(0., 1., len)

    ils_spl = Spline1D(ils.Loc[:,2], R2r(ils.Loc[:,1]))
    rILS = ils_spl(μcol)
    yILS = rILS .* μcol
    xILS = sqrt(rILS.^2 - yILS.^2)

    rIRS = 1.+sqrt(1-crd.a^2 * μcol.^2)
    yIRS = rIRS .* μcol
    xIRS = sqrt(rIRS.^2 - yIRS.^2)

    rhz = rmin
    yhz = rhz * μcol
    xhz = sqrt(rhz^2 - yhz.^2)

    levels = linspace(0.005, Umax, cnum)
    figure(figsize=(5, 6))
    contour(Uxy, levels, extent = (0, xmax, 0, ymax), colors = "k")
    plot(xIRS, yIRS, "r-")
    plot(xILS, yILS,  "k--")
    fill_between(xhz, 0., yhz, color = "black")
    xlabel(L"$X/M$", fontsize = 20)
    ylabel(L"$Z/M$", fontsize = 20)
    tight_layout()
    #savefig("f1.pdf")
end

function Fsq(U::Array{Float64, 2}, crd::Cord, grd::Grid, Ω_I::Ω_and_I)
    idx = crd.idx_r2
    ∂1U = zeros(idx);  ∂1U[2:idx] = (U[1, 2:idx] - U[1, 1:idx-1]) ./ crd.δR
    ∂2U = (U[2, 1:idx] - U[1, 1:idx]) ./ crd.δμ

    Rcol = crd.R[1, 1:idx]
    ∂rU  = ∂1U .* (1-Rcol).^2; ∂μU = ∂2U

    μ = 0.
    r = crd.rcol[1:idx]
    Δ = r.^2 - 2r + crd.a^2
    Σ = r.^2 + crd.a^2 * μ.^2
    β = Δ .* Σ + 2r .*(r.^2+crd.a^2)

    Ωspl = Ω_I.Ωspl
    Ucol = U[1,1:idx]
    Icol = 2*Ωspl(Ucol).*Ucol
    κcol = grd.κ[1,1:idx]

    B2mE2 =                -κcol .* (Δ .* ∂rU.^2 + ∂μU.^2)./Σ + Icol.^2
    B2pE2 = (κcol  + Δ .* Σ ./β ).* (Δ .* ∂rU.^2 + ∂μU.^2)./Σ + Icol.^2
    fsq   = B2mE2./B2pE2

    fsq_exp   = 0.*exp(-(r-crd.rmin)/0.1)
    B2mE2_exp = B2pE2.*fsq_exp
    δ_B2mE2   = B2mE2 - B2mE2_exp

    fsq2_spl = Spline1D(Ucol[2:end], (fsq[2:end]-fsq_exp[2:end]).^2)
    fsq2_avg = integrate(fsq2_spl, Ucol[2], Ucol[end])/(Ucol[end]-Ucol[2])

    # plot(Ucol, Icol.^2, "k")
    # plot(Ucol, κcol .* (Δ .* ∂rU.^2 + ∂μU.^2)./Σ, "r")
    # plot(Ucol, κcol .* (Δ .* ∂rU.^2)./Σ, "b--")
    # plot(Ucol, κcol .* (∂μU.^2)./Σ, "r--")
    plot(Ucol, B2mE2)
    plot(Ucol, B2mE2_exp, "--")
    return Ucol, δ_B2mE2, fsq, fsq2_avg
end

function Fsq_only(U::Array{Float64, 2}, crd::Cord, grd::Grid, Ω_I::Ω_and_I)
    idx = crd.idx_r2
    ∂1U = zeros(idx);  ∂1U[2:idx] = (U[1, 2:idx] - U[1, 1:idx-1]) ./ crd.δR
    ∂2U = (U[2, 1:idx] - U[1, 1:idx]) ./ crd.δμ

    Rcol = crd.R[1, 1:idx]
    ∂rU  = ∂1U .* (1-Rcol).^2; ∂μU = ∂2U

    μ = 0.
    r = crd.rcol[1:idx]
    Δ = r.^2 - 2r + crd.a^2
    Σ = r.^2 + crd.a^2 * μ.^2
    β = Δ .* Σ + 2r .*(r.^2+crd.a^2)

    Ωspl = Ω_I.Ωspl
    Ucol = U[1,1:idx]
    Icol = 2*Ωspl(Ucol).*Ucol
    κcol = grd.κ[1,1:idx]

    B2mE2 = -κcol .* (Δ .* ∂rU.^2 + ∂μU.^2)./Σ + Icol.^2
    B2pE2 = (κcol  + Δ .* Σ ./β ).* (Δ .* ∂rU.^2 + ∂μU.^2)./Σ + Icol.^2

    fsq_exp   = 0.*exp(-(r-crd.rmin)/0.1)
    B2mE2_exp = B2pE2.*fsq_exp

    δ_B2mE2 = B2mE2 - B2mE2_exp
    return Ucol, δ_B2mE2
end


function Edot(Ubm, Ωbm)
    Ibm = 2*Ωbm.*Ubm
    xspl = Spline1D(Ubm, Ωbm.*Ibm)
    return integrate(xspl, Ubm[1], Ubm[end]) * pi
end
