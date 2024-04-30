####################################################################################
## Estimate Coeffients in Aggregate Capital Transition in Krusell & Smith (1998)  ##
####################################################################################

function KS_OLS(m, kbar_sim, z_sim_index)
    """
    === Estimate {a0, a1, b0, b1} ===
    <input>
    ・m:
    ・kbar_sim (t = T0+1 - T_sim)
    ・z_sim_index (t = 1 - T_sim) 
    """
    T0 = m.T0
    
    ind_1 = findall(z_sim_index[T0+1:end-1] .== 1.0)
    ind_2 = findall(z_sim_index[T0+1:end-1] .== 2.0)
    
    X1 = [ones(length(ind_1)) log.(kbar_sim[ind_1])]
    X2 = [ones(length(ind_2)) log.(kbar_sim[ind_2])]

    Y1 = log.(kbar_sim[ind_1 .+ 1]);
    Y2 = log.(kbar_sim[ind_2 .+ 1]);

    a0_hat, a1_hat = (X2'*X2)^(-1) * (X2'Y2)
    b0_hat, b1_hat = (X1'*X1)^(-1) * (X1'Y1)

    return a0_hat, a1_hat, b0_hat, b1_hat 
end


function KS_R2(m,kbar_sim,z_sim_index,a0,a1,b0,b1)
    """
    === Compute R^2 ===
    """

    T0 = m.T0
    
    a = [a0; a1];
    b = [b0; b1];

    ind_1 = findall(z_sim_index[T0+1:end-1] .== 1.0)
    ind_2 = findall(z_sim_index[T0+1:end-1] .== 2.0)
    
    X1 = [ones(length(ind_1)) log.(kbar_sim[ind_1])]
    X2 = [ones(length(ind_2)) log.(kbar_sim[ind_2])]

    Y1 = log.(kbar_sim[ind_1 .+ 1]);
    Y2 = log.(kbar_sim[ind_2 .+ 1]);

    e1 = Y1 .- X1 * b;
    e2 = Y2 .- X2 * a;
    sig_bad = (e1'*e1) / length(ind_1);
    sig_good = (e2'*e2) / length(ind_2);

    i_vec1 = ones(length(ind_1));
    i_vec2 = ones(length(ind_2));

    I1 = Matrix(I,length(ind_1),length(ind_1));
    I2 = Matrix(I,length(ind_2),length(ind_2));

    M1 = I1 - i_vec1 * (i_vec1' * i_vec1)^(-1) * i_vec1'
    M2 = I2 - i_vec2 * (i_vec2' * i_vec2)^(-1) * i_vec2'

    R2_good = (a' * X2' * M2 * X2 * a) / (Y2' * M2 * Y2)
    R2_bad =  (b' * X1' * M1 * X1 * b) / (Y1' * M1 * Y1)

    return R2_good, R2_bad, sig_good, sig_bad
end
