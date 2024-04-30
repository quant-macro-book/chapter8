################################################
## HH Optimization in Krusell & Smith (1998)  ##
################################################

include("KS_markov.jl")

struct Model{TI<:Integer,TF<:AbstractFloat}
    # primitive params
    alpha::TF
    delta::TF
    beta::TF
    lbar::TF
    
    # shock grid
    Ne::TI
    egrid::Array{TF,1} # idiosyncratic shock
    Nz::TI
    zgrid::Array{TF,1} # aggregate shock: {z_{b}, z_{g}}
    Nu::TI
    ugrid::Array{TF,1} # unemployment rate depending on the aggregate shock: {u_{b}, u_{g}}
    pi_z::Array{TF,2} # Π(z'|z)
    pi_e::Array{TF,2} # Π(e'|z',z,e)
    
    # capital grid
    Nm::TI
    mmin::TF
    mmax::TF
    mgrid::Array{TF,1} # approximate aggregation

    Nk::TI
    Nk_dist::TI
    kmin::TF
    kmax::TF
    kgrid::Array{TF,1}
    kgrid_dist::Array{TF,1}

    # parameters for Simulation
    T_sim::TI
    T0::TI
end

function Construct(
    # primitive params
    alpha = 0.36,
    delta = 0.025,
    beta = 0.99,
    lbar = 0.3271,

    # shock grid
    Ne = 2,
    egrid = [0.0001, 1.0],
    Nz = 2,
    zgrid = [0.99,1.01],
    Nu = 2,
    ugrid = [0.1, 0.04],

    # capital grid
    Nm = 5,
    mmin = 9.0,
    mmax = 14.0,
    Nk = 35,
    Nk_dist = 1000,
    kmin = 0.0,
    kmax = 75.0,

    # parameters for Simulation
    T_sim = 6000,
    T0 = 1000
    )
    
    # shock grid
    pi_z, pi_e = KS_markov(ugrid, Ne, Nz);

    # capital grid
    mgrid = collect(range(mmin, mmax, length= Nm));
    kgrid = exp.(range(0,log(kmax-kmin+1),length=Nk)) .+ (-1+kmin) # 1st order log spaced grid
    kgrid_dist = collect(range(kmin,kmax,length=Nk_dist))

    
    
    
    return Model(alpha, delta, beta, lbar, Ne, egrid, Nz, zgrid, Nu, ugrid, pi_z, pi_e, Nm, mmin, mmax, mgrid, Nk, Nk_dist, kmin, kmax, kgrid, kgrid_dist, T_sim, T0)
end


# VFI 
@inline function KS_vfi(m::Model, a0::Float64, a1::Float64, b0::Float64, b1::Float64, v0::Array{Float64,4}, tol=1e-7, maxit=5000, diff=1, num=0)
    """
    === Solving Household Problem by VFI in Krusell & Smith (1998) ===
    <input>
    ・m:
    ・a0:
    ・a1:
    ・b0: 
    ・b1:
    ・v0:
    ・tol:
    ・maxit:
    ・diff:
    ・num:
    <output>
    ・Tv:
    ・g:
    ・c: 
    """

    alpha, delta, beta, lbar, Ne, egrid, Nz, zgrid, Nu, ugrid = m.alpha, m.delta, m.beta, m.lbar, m.Ne, m.egrid, m.Nz, m.zgrid, m.Nu, m.ugrid
    pi_z, pi_e, Nm, mgrid, Nk, kmin, kmax, kgrid = m.pi_z, m.pi_e, m.Nm, m.mgrid, m.Nk, m.kmin, m.kmax, m.kgrid

    # (1) Initialization
    v = v0
    g = similar(v0)
    c = similar(v0)
    Tv = similar(v0)
    m_update = [b0 b1; a0 a1]
    
    # (2) VFI
    while (diff > tol) & (num < maxit)

        Threads.@threads for (im, mnow) in collect(enumerate(mgrid))
            Threads.@threads for (iz, znow) in collect(enumerate(zgrid))

                # perice 
                r = alpha * znow * (mnow^(alpha-1.0)) * (((1.0-ugrid[iz])*lbar)^(1.0-alpha)) - delta
                w = (1.0-alpha) * znow * (mnow^alpha) * (((1.0-ugrid[iz])*lbar)^(-alpha))

                Threads.@threads for (ik, know) in collect(enumerate(kgrid))
                    Threads.@threads for (ie, enow) in collect(enumerate(egrid))

                        # update kbar
                        mprime = exp(m_update[iz,1] + m_update[iz,2]*log(mnow)) # kbar'

                        # expectation and interpolation
                        ev = EV(m,ie,iz,v)
                        v_spline = Spline2D(kgrid,mgrid,ev,kx=3,ky=3,s=0.0)

                        # solve bellman equation
                        income = (1+r)*know + w*enow*lbar
                        ub = min(income, kmax)

                        res = optimize(x -> Bellman(x,mprime,income,v_spline),kmin,ub,GoldenSection())
                        g[ik,ie,iz,im], Tv[ik,ie,iz,im] = res.minimizer, res.minimum
                        Tv[ik,ie,iz,im] = -Tv[ik,ie,iz,im]
                        c[ik,ie,iz,im] = income - g[ik,ie,iz,im]

                    end
                end
            end
        end

        # check convergence
        diff = maximum(abs.(Tv-v))
        num += 1 
        # @printf "iteration =%4d, diff=%8.6f \n" num diff 
        # flush(stdout)
        v = copy(Tv) # update
    end

    return Tv, g, c
end

@inline function EV(m::Model, ie::Int64, iz::Int64, v::Array{Float64,4})
    """
    === Expectation of Value Function Next Period ===
    """
    Ne, Nz, Nm, Nk, pi_z, pi_e = m.Ne, m.Nz, m.Nm, m.Nk, m.pi_z, m.pi_e

    ev = zeros(Nk, Nm)

    @inbounds for izz in 1:Nz

        temp1 = pi_z[iz,izz]
        
        @inbounds for iee in 1:Ne
            
            temp2 = pi_e[2*(iz-1)+ie, 2*(izz-1)+iee] / sum(pi_e[2*(iz-1)+ie, 2*(izz-1)+1:2*(izz-1)+2])
            ev += temp1 * temp2 * v[:,iee,izz,:]
            
        end
    end

    return ev
end


@inline function Bellman(kprime,mprime,income,v_spline)
    """
    === RHS of Bellman Equation ===
    """
    beta = m.beta

    return -log(income-kprime) - beta * v_spline(kprime,mprime)

end