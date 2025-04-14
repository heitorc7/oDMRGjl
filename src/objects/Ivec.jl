function make_ivec(L::Int, sites)
    # Initialize MPS with the given sites
    ψ = MPS(sites)
    
    # First site
    A = ψ[1]
    inds_i = inds(A)
    B = ITensor(inds_i[1], inds_i[2])
    B[inds_i[1](1), inds_i[2](1)] = 1
    B[inds_i[1](2), inds_i[2](1)] = 0
    B[inds_i[1](3), inds_i[2](1)] = 0
    B[inds_i[1](4), inds_i[2](1)] = 1
    ψ[1] = B
    
    # Bulk sites
    for i in 2:L-1
        A = ψ[i]
        inds_i = inds(A)
        B = ITensor(inds_i[1], inds_i[2], inds_i[3])
        
        B[inds_i[1](1), inds_i[2](1), inds_i[3](1)] = 1
        B[inds_i[1](1), inds_i[2](2), inds_i[3](1)] = 0
        B[inds_i[1](1), inds_i[2](3), inds_i[3](1)] = 0
        B[inds_i[1](1), inds_i[2](4), inds_i[3](1)] = 1
        
        ψ[i] = B
    end
    
    # Last site
    A = ψ[L]
    inds_i = inds(A)
    B = ITensor(inds_i[1], inds_i[2])
    B[inds_i[1](1), inds_i[2](1)] = 1
    B[inds_i[1](1), inds_i[2](2)] = 0
    B[inds_i[1](1), inds_i[2](3)] = 0
    B[inds_i[1](1), inds_i[2](4)] = 1
    ψ[L] = B
    
    return ψ
end