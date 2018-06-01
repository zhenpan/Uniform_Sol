using Dierckx

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
    Res  = zeros(μlen, Rlen)
    dU   = zeros(U)
    Upls = zeros(crd.μlen)
    Umns = zeros(crd.μlen)

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

        for j = 2:μlen-1                                 #odd run for even j, even run for odd j
            for l = 3-mod(j,2):2:idx_xbd[j]-1
                if lsn_map[j,l] == 0
                    Res[j,l] = a[j,l]*U[j+1,l] + b[j,l]*U[j-1,l] + c[j,l]*U[j,l+1] + d[j,l]*U[j,l-1] + ee[j,l]*U[j,l] - f[j,l]
                    dU[j,l]  = -omega*Res[j,l]/ee_rgl[j,l]
                end
            end
        end

        U += dU
        U = USmooth!(U, lsn, crd)       #lsn bounds, updated via interior points
    end
    return U, Res, dU
end


function Bounds!(U::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I, lsn::LS_neighbors, bc_eqt::BC_eqt)
    idx_r2  = crd.idx_r2
    idx_bd  = crd.idx_xbd[1]

    #horizon and inf r boundary values
    U[:,1]   = U[:,2]             # computation friendly BC on horizon
    U[:,end] = U[:,end-1]         # inf r boundary is in fact not used due to xbd

    Utmp = U[2, 1:idx_bd] - bc_eqt.∂μU*crd.δμ                   # equator boundary
    U[1, 1:idx_bd] = Utmp;  U[1,idx_r2+1] = 1.002*U[1,idx_r2]

    δR  = crd.δR
    U_H = U[1, idx_r2]*(crd.Rcol[idx_r2+1]-r2R(2.0))/δR + U[1, idx_r2+1]*(r2R(2.0)-crd.Rcol[idx_r2])/δR
    return U, U_H
end

function BC_gen(U::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I; BC_opt = 0, Isf = 5.)
    idx_r2  = crd.idx_r2
    idx_bd  = crd.idx_xbd[1]
    rmin    = crd.rmin
    ∂μU     = zeros(idx_bd)

    if BC_opt==0
        ∂μU[1:idx_r2] = -1.*(crd.rcol[1:idx_r2]/rmin).^2.8
    else
        ∂μU[1:idx_r2] = (U[2,1:idx_r2]-U[1,1:idx_r2])/crd.δμ

        Ucol= U[1,1:idx_r2]
        rcol= crd.rcol[1:idx_r2]
        Ubm = linspace(0., Ucol[end], 128)
        Ibm = 2*Ω_I.Ωspl(Ubm).*Ubm

        Ispl = Spline1D(Ubm,Ibm)
        iip  = Ispl(Ucol).*derivative(Ispl, Ucol)
        IIp  = Ω_I.IIpspl(Ucol)

        ∂μUnew  = ∂μU[1:idx_r2] + Isf*(IIp-iip)
        pmodel(x, p) = ( p[1] + p[2] .* x + p[3] .* x.^2 + p[4] .* x.^3 + p[5] .* x.^4 + p[6] .* x.^5 )
        pfit   = curve_fit(pmodel, rcol, ∂μUnew, [0., 0., 0., 0., 0., 0.])
        ∂μU[1:idx_r2] = pmodel(rcol, pfit.param)
    end

    ∂μU[idx_r2+1:idx_bd] = ∂μU[idx_r2]*exp(-(crd.rcol[idx_r2+1:idx_bd]-2.).^2/(2*0.05^2))
    return BC_eqt(∂μU)
end

function Init(crd::Cord, mtr::Geom; xbd = 4.0)
        z = crd.r .* crd.μ
        x = sqrt(crd.r.^2 - z.^2)
        U = x.^2 + 0.65acos(crd.μ).*exp(-(crd.r-0.8).^2).*exp(-z)
        U_H = U[1, crd.idx_r2]

        #initialize Ω_and_I
        Ubm = collect(linspace(0., U_H, 2048)); Ωbm = zeros(Ubm)
        for i = 1:length(Ubm)
            #Ωbm[i] = (Ubm[i] < U_H) ? 0.5*crd.Ω_H*(cos(pi/2*Ubm[i]/U_H).^2) : 0.
            Ωbm[i] = (Ubm[i] < U_H) ? 0.5*crd.Ω_H*(1-Ubm[i]/U_H) : 0.
        end

        Ibm  = 2*Ωbm.*Ubm
        Ωspl = Spline1D(Ubm, Ωbm, bc = "zero")
        Ispl = Spline1D(Ubm, Ibm, bc = "zero")
        Ipbm = derivative(Ispl, Ubm)
        IIpspl= Spline1D(Ubm, Ibm.*Ipbm, bc = "zero")
        Ω_I   = Ω_and_I(U, crd, Ωspl, Ispl, IIpspl)

        return U, Ω_I, U_H
end

#
# function BC_init(U::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I)
#     #equator bounds ( within and beyond r2)
#     idx_r2  = crd.idx_r2
#     idx_bd  = crd.idx_xbd[1]
#
#     # Ispl = Ω_I.Ispl     #I_solver(Ω_I)
#     # Ωspl = Ω_I.Ωspl
#     # Uhe  = 2.5 #U[1,1]       #U in the horizon/equator cornor
#     # Ωhe  = Ωspl(Uhe)
#     # Ihe  = Ispl(Uhe)
#
#     Ω_H  = crd.Ω_H
#     rmin = crd.rmin
#
#     ∂μ         = zeros(idx_bd)
#     ∂μ[1]      = -1. #0.5*rmin* Ihe/(Ωhe-Ω_H)                     # obtain from Znajek Condition, ∂μ shoule be negative
#     ∂μ[1:idx_r2] = ∂μ[1]*(crd.rcol[1:idx_r2]/rmin).^2.8 - 1.0*exp(-(crd.rcol[1:idx_r2]-1.85).^2/(2*0.2^2))-0.2*exp(-(crd.rcol[1:idx_r2]-1.5).^2/(2*0.15^2)) +0.2*exp(-(crd.rcol[1:idx_r2]-rmin).^2/(2*0.5^2)) # initial guess
#     ∂μ[idx_r2+1:idx_bd] = ∂μ[idx_r2]*exp(-(crd.rcol[idx_r2+1:idx_bd]-2.).^2/(2*0.05^2))
#
#
#     # U[1, idx_r2+1:idx_bd]= U[2, idx_r2+1:idx_bd] - ∂μ[idx_r2+1:idx_bd]*crd.δμ
#     # U[1, idx_r2]         = U[1, idx_r2+1]
#     # U[1,1]               = U[2,1]-∂μ[1]*crd.δμ
#     #
#     # A = (U[1,1]-U[1,idx_r2])/(rmin-2)^2
#     # U[1,2:idx_r2-1] = A*(crd.rcol[2:idx_r2-1]-2).^2 + U[1, idx_r2]
#     return ∂μ
# end
