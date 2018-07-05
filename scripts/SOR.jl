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


function Bounds!(U::Array{Float64,2}, crd::Cord, Ω_I::Ω_and_I, lsn::LS_neighbors, bc_eqt::BC_eqt)
    idx_r2  = crd.idx_r2
    idx_bd  = crd.idx_xbd[1]

    #horizon and inf r boundary values
    U[2:end,1] = U[2:end,2]         # computation friendly BC on horizon
    U[:,end]   = U[:,end-1]         # inf r boundary is in fact not used due to xbd


    U[1, 1:idx_r2] = bc_eqt.Ueqt
    U[1, idx_r2+1] = U[1,idx_r2]
    U[1, idx_r2+2] = U[1,idx_r2]
    U[1, idx_r2+3:idx_bd]   = U[2, idx_r2+3:idx_bd]               # equator boundary

    U = USmooth!(U, lsn, crd)

    δR  = crd.δR
    U_H = U[1, idx_r2]*(crd.Rcol[idx_r2+1]-r2R(2.0))/δR + U[1, idx_r2+1]*(r2R(2.0)-crd.Rcol[idx_r2])/δR
    return U, U_H
end

function BC_gen(U::Array{Float64,2}, crd::Cord, grd::Grid, Ω_I::Ω_and_I; BC_opt = 0, Isf = 5., Isf_BE = 0.)
    idx_r2  = crd.idx_r2
    rmin    = crd.rmin
    Ueqt    = zeros(idx_r2)

    if BC_opt==0
        ∂μU[1:idx_r2] = -1.5*(crd.rcol[1:idx_r2]/rmin).^2.8 #-1.*(crd.rcol[1:idx_r2]/rmin).^2.8
    else
        ∂μU = (U[2,1:idx_r2]-U[1,1:idx_r2])/crd.δμ
        Ucol= U[1,1:idx_r2]
        rcol= crd.rcol[1:idx_r2]
        Ubm = linspace(0., Ucol[end], 128)
        Ibm = 2*Ω_I.Ωspl(Ubm).*Ubm

        Ispl = Spline1D(Ubm,Ibm)
        iip  = Ispl(Ucol).*derivative(Ispl, Ucol)
        IIp  = Ω_I.IIpspl(Ucol)

        Ueqt, B2mE2 = Fsq_only(U, crd, grd, Ω_I)

        ∂μUnew  = ∂μU[1:idx_r2] + Isf*(IIp-iip) - Isf_BE*(1.-rcol/rcol[end]).*B2mE2
        pmodel(x, p) = ( p[1] + p[2] .* x + p[3] .* x.^2 + p[4] .* x.^3 + p[5] .* x.^4 + p[6] .* x.^5 )
        pfit   = curve_fit(pmodel, rcol, ∂μUnew, [0., 0., 0., 0., 0., 0.])
        ∂μU[1:idx_r2] = pmodel(rcol, pfit.param)

        (xmax, idx_max) = findmax(∂μU[1:idx_r2])
        ∂μU[1:idx_max]  = xmax
    end

    return BC_eqt(Ueqt)
end

function Init(crd::Cord, mtr::Geom, Ω_par::Array{Float64}; xbd = 4.0)
        z = crd.r .* crd.μ
        x = sqrt(crd.r.^2 - z.^2)
        U = x.^2

        #initialize equator values
        rmin = crd.rmin
        xeq  = vcat( linspace(0., 4., 512), logspace(log10(4.01), log10(crd.rcol[end]), 512) )
        Ueq  = zeros(xeq); U_he = 1.6; U_H = 4.4

        A = (U_H-U_he + 0.8*(2-rmin)^2)/(2-rmin)

        for i = 1:length(Ueq)
            if xeq[i] < rmin
                Ueq[i] =  U_he* (1. - (xeq[i]/rmin-1)^2)
            elseif xeq[i] < 2.
                Ueq[i] = A*(xeq[i]-rmin)- 0.8*(xeq[i]-rmin)^2 +U_he
            elseif xeq[i] < xbd^2
                Ueq[i] = xeq[i]^2 + (U_H-4.)*exp(2^2-xeq[i]^2)
            else
                Ueq[i] = xeq[i].^2
            end
        end
        Uspl = Spline1D(xeq, Ueq, k=1)

        #initialize all grid points via interpolation
        for j = 1:crd.μlen
            for l = 1:crd.Rlen
                U[j,l] = U[j,l] + (Uspl(x[j,l])-U[j,l])*exp(-3*z[j,l]) #.*exp(-crd.μ[j,l]^3)
            end
        end

        U_H = U[1, crd.idx_r2]

        #initialize Ω_and_I
        Ubm = collect(linspace(0., U_H, 2048))
        Ωbm = Ω_fnc(crd.Ω_H, Ω_par, Ubm/U_H)

        Ibm  = 2*Ωbm.*Ubm
        Ωspl = Spline1D(Ubm, Ωbm, bc = "zero")
        Ispl = Spline1D(Ubm, Ibm, bc = "zero")
        Ipbm = derivative(Ispl, Ubm)
        IIpspl= Spline1D(Ubm, Ibm.*Ipbm, bc = "zero")
        Ω_I   = Ω_and_I(U, crd, Ωspl, Ispl, IIpspl)

        return U, Ω_I, U_H
end
