#################################################################
## Simulation for Aggregate Capital in Krusell & Smith (1998)  ##
#################################################################

function sim_kbar(m::Model, g::Array{Float64,4}, z_sim_index::Array{Float64,1}, mu0::Array{Float64,2})
    """
    === Simulate \bar{k} ===
    <input>
    ・m
    ・g
    ・z_sim: simulated aggregate shock array
    ・mu0: 
    """

    kgrid, kgrid_dist, Nk, Nk_dist, Ne, zgrid, Nz, mgrid, Nm = m.kgrid, m.kgrid_dist, m.Nk, m.Nk_dist, m.Ne, m.zgrid, m.Nz, m.mgrid, m.Nm
    pi_z, pi_e, T_sim, T0 = m.pi_z, m.pi_e, m.T_sim, m.T0

    kbar_sim = zeros(T_sim)
    mu = zeros(T_sim, Nk_dist, Ne)
    mu[1,:,:] .= copy(mu0)
    kbar_sim[1] = sum(kgrid_dist .* sum(mu[1,:,:], dims=2)) # \bar{k_{t}}

    for t in 1:T_sim-1

        iz = Int(z_sim_index[t]) # aggregate shock this period
        izz = Int(z_sim_index[t+1]) # aggregate shock next period
        
        # interpolation (1): \bar{k}
        g_dist_pre = zeros(Nk, Ne) 
        g_dist = zeros(Nk_dist, Ne)

        Threads.@threads for ie in 1:Ne
            Threads.@threads for ik in 1:Nk
                
                g_interp_pre = Spline1D(mgrid, g[ik,ie,iz,:], k=1, bc="extrapolate")
                g_dist_pre[ik,ie] = g_interp_pre(kbar_sim[t])

            end
        end

        # interpolation (2): kgrid_dist
        Threads.@threads for ie in 1:Ne

            g_interp = Interpolator(kgrid, g_dist_pre[:,ie])
            g_dist[:,ie] = g_interp.(kgrid_dist)

        end

        @inbounds for ik in 1:Nk_dist # asset this period
            @inbounds for ie in 1:Ne  # idiosyncratic shock this period

                index_k = max(searchsortedlast(kgrid_dist,g_dist[ik,ie]),1)
                wj = (kgrid_dist[index_k+1] - g_dist[ik, ie]) /  (kgrid_dist[index_k+1] - kgrid_dist[index_k])
                
                # update the distribution
                @inbounds for iee in 1:Ne # ideosyncratic shock next period
                    
                    mu[t+1, index_k, iee] += wj * (pi_e[2*(iz-1)+ie, 2*(izz-1)+iee] / sum(pi_e[2*(iz-1)+ie,2*(izz-1)+1:2*(izz-1)+2])) * mu[t, ik, ie] 
                    mu[t+1, index_k+1, iee] += (1-wj) * (pi_e[2*(iz-1)+ie, 2*(izz-1)+iee] / sum(pi_e[2*(iz-1)+ie,2*(izz-1)+1:2*(izz-1)+2])) * mu[t, ik, ie]

                end
            end
        end

        #mu[t+1,:,:] /= sum(mu[t+1,:,:])

        kbar_sim[t+1] = sum(kgrid_dist .* sum(mu[t+1,:,:], dims=2)) # \bar{k_{t}}

    end

    return kbar_sim[T0+1:end]
end


function KS_stationary(g,m::Model,mu,tol=1e-8,maxit=5000,distance=1.0,num=1)
    """
    === Calculating Stationary Distribution by Aiyagari Approach === 
    """
    kgrid, kgrid_dist, Nk_dist = m.kgrid, m.kgrid_dist, m.Nk_dist
    Nm, Ne, egrid, pi_e = m.Nm, m.Ne, m.egrid, m.pi_e
    
    #(1) Re-generate a new grid, and interpolation
    g_dist = zeros(Nk_dist, Ne)

    @inbounds for ie in 1:Ne
        #g_interp = LinearInterpolation(agrid, g[:, is])
        g_interp = Interpolator(kgrid, g[:,ie,2,end])
        g_dist[:, ie] .= g_interp.(kgrid_dist)
    end

    P = m.pi_e[3:4,3:4] ./ sum(m.pi_e[3:4,3:4], dims=2) # 今期と来期の aggregate state が "good" を仮定

    #(2) Solving stationary distribution
    while (num < maxit) & (distance > tol)
        Tmu = zeros(Nk_dist, Ne)
        @inbounds for ik in 1:Nk_dist # asset this period
            @inbounds for ie in 1:Ne  # state this period
                index_k = max(searchsortedlast(kgrid_dist,g_dist[ik,ie]),1)
                wj = (kgrid_dist[index_k+1] - g_dist[ik,ie]) /  (kgrid_dist[index_k+1] - kgrid_dist[index_k])
                # update the distribution
                @inbounds for iee in 1:Ne #state next period
                    Tmu[index_k, iee] +=  wj * P[ie,iee] * mu[ik, ie]
                    Tmu[index_k+1, iee] += (1-wj) * P[ie,iee] * mu[ik,ie]
                end
            end
        end

        #check the convergence
        distance = maximum(abs.(Tmu-mu))
        num += 1 
        #@printf "iteration =%4d, diff=%8.6f \n" num distance 
        mu = copy(Tmu) #update
        if num == maxit
            error("distribution does not converge.")
        end
    end
    return mu
end