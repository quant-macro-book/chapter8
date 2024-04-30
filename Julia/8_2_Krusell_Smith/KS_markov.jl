#################################################################
## Markov Chain of Idiosyncratic & Aggregate Shock in KS(1998) ##
## Π(z'|z) & Π(e'|z',z,e)                                      ##
#################################################################

function KS_markov(ugrid, Ne, Nz)
    """
    ==============================
    Create Π(z'|z) & Π(e'|z',z,e)
    ==============================
    <input>
    ・ugrid:
    ・Ne: 
    ・Nz: 
    <output>
    ・pi_e: Π(e'|z',z,e)
    ・pi_z: Π(z'|z)
    """

    # ugrid, Ne, Nz =  m.ugrid, m.Ne, m.Nz

    pi_e = zeros(Ne*Nz, Ne*Nz)
    pi_z = zeros(Nz, Nz)

    # (1) create Π(z'|z)
    Az = [8 0 0 0; 0 8 0 0; 1 0 0 1; 0 1 1 0]
    bz = [7; 7; 1; 1]
    xz = Az^-1 * bz
    
    pi_z[1,1] = xz[1]
    pi_z[2,2] = xz[2]
    pi_z[1,2] = xz[3]
    pi_z[2,1] = xz[4]

    # (2) create Π(e'|z',z,e)
    Ap = zeros((Ne*Nz)^2,(Ne*Nz)^2)

    Ap[1,11] = 1.5
    Ap[2,1] = 2.5
    Ap[3,1] = -1.25*xz[4]
    Ap[3,9] = xz[2]
    Ap[4,3] = xz[1]
    Ap[4,11] = -0.75*xz[3]

    Ap[5,11] = 1.0
    Ap[5,12] = 1.0
    Ap[6,15] = 1.0
    Ap[6,16] = 1.0
    Ap[7,9] = 1.0
    Ap[7,10] = 1.0
    Ap[8,13] = 1.0
    Ap[8,14] = 1.0
    Ap[9,3] = 1.0
    Ap[9,4] = 1.0
    Ap[10,7] = 1.0
    Ap[10,8] = 1.0
    Ap[11,1] = 1.0
    Ap[11,2] = 1.0
    Ap[12,5] = 1.0
    Ap[12,6] = 1.0

    Ap[13,11] = ugrid[2]
    Ap[13,15] = 1-ugrid[2]
    Ap[14,9] = ugrid[2]
    Ap[14,13] = 1-ugrid[2]
    Ap[15,3] = ugrid[1]
    Ap[15,7] = 1-ugrid[1]
    Ap[16,1] = ugrid[1]
    Ap[16,5] = 1-ugrid[1]

    bp = [0.5*xz[1]; 1.5*xz[2]; 0.0; 0.0; xz[1]; xz[1]; xz[4]; xz[4]; xz[3]; xz[3]; xz[2]; xz[2]; ugrid[2]*xz[1]; ugrid[1]*xz[4]; ugrid[2]*xz[3]; ugrid[1]*xz[2]]
    xp = Ap^-1 * bp

    for i in 1:Ne*Nz # 行
        for j in 1:Ne*Nz # 列
            pi_e[i,j] = xp[(Nz*Ne)*(i-1) + j]
        end
    end

    return pi_z, pi_e
end


function markov_sim(m)
    """
    === Simulation for Aggregate Shock ===
    """

    zgrid, pi_z, T_sim = m.zgrid, m.pi_z, m.T_sim


    # (1) preparation
    e = rand(Uniform(0,1),T_sim) # (T×1) vector of random draws of uniform distribution
    z_sim = zeros(T_sim) # array for simulated value of z (0.99 or 1.01)
    z_sim_index = similar(z_sim) # array for simulated "index" value of z (1 or 2)
    z_sim_index[1] = 2 # index for the initial value of z
    z_sim[1] = zgrid[Int(z_sim_index[1])] # initial value is bad state

    pi_z_hat = cumsum(pi_z, dims=2)

    # (2) simulations
    for t in 1:T_sim-1
        pi_z_hat_now = pi_z_hat[Int(z_sim_index[t]),:]
        z_sim_index[t+1] = sum(e[t] .- pi_z_hat_now .>= 0)
        z_sim_index[t+1] = min(Int(z_sim_index[t+1]+1), 2)
    
        z_sim[t+1] = zgrid[Int(z_sim_index[t+1])]
    end

    return z_sim_index, z_sim
end